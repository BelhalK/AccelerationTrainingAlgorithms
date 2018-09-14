from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import six
import tensorflow as tf

from edward.inferences.variational_inference import VariationalInference
from edward.models import RandomVariable
from edward.util import copy, get_descendants

try:
  from edward.models import Normal
  from tensorflow.contrib.distributions import kl_divergence
except Exception as e:
  raise ImportError("{0}. Your TensorFlow version is not supported.".format(e))


class test(VariationalInference):
  """Variational inference with the KL divergence
  $\\text{KL}( q(z; \lambda) \| p(z \mid x) ).$
  This class minimizes the objective by automatically selecting from a
  variety of black box inference techniques.
  #### Notes
  `KLqp` also optimizes any model parameters $p(z \mid x;
  \\theta)$. It does this by variational EM, maximizing
  $\mathbb{E}_{q(z; \lambda)} [ \log p(x, z; \\theta) ]$
  with respect to $\\theta$.
  In conditional inference, we infer $z$ in $p(z, \\beta
  \mid x)$ while fixing inference over $\\beta$ using another
  distribution $q(\\beta)$. During gradient calculation, instead
  of using the model's density
  $\log p(x, z^{(s)}), z^{(s)} \sim q(z; \lambda),$
  for each sample $s=1,\ldots,S$, `KLqp` uses
  $\log p(x, z^{(s)}, \\beta^{(s)}),$
  where $z^{(s)} \sim q(z; \lambda)$ and $\\beta^{(s)}
  \sim q(\\beta)$.
  The objective function also adds to itself a summation over all
  tensors in the `REGULARIZATION_LOSSES` collection.
  """
  def __init__(self, latent_vars=None, data=None):
    """Create an inference algorithm.
    Args:
      latent_vars: list of RandomVariable or
                   dict of RandomVariable to RandomVariable.
        Collection of random variables to perform inference on. If
        list, each random variable will be implictly optimized using a
        `Normal` random variable that is defined internally with a
        free parameter per location and scale and is initialized using
        standard normal draws. The random variables to approximate
        must be continuous.
    """
    if isinstance(latent_vars, list):
      with tf.variable_scope(None, default_name="posterior"):
        latent_vars_dict = {}
        continuous = \
            ('01', 'nonnegative', 'simplex', 'real', 'multivariate_real')
        for z in latent_vars:
          if not hasattr(z, 'support') or z.support not in continuous:
            raise AttributeError(
                "Random variable {} is not continuous or a random "
                "variable with supported continuous support.".format(z))
          batch_event_shape = z.batch_shape.concatenate(z.event_shape)
          loc = tf.Variable(tf.random_normal(batch_event_shape))
          scale = tf.nn.softplus(
              tf.Variable(tf.random_normal(batch_event_shape)))
          latent_vars_dict[z] = Normal(loc=loc, scale=scale)
        latent_vars = latent_vars_dict
        del latent_vars_dict

    super(test, self).__init__(latent_vars, data)

  def initialize(self, n_samples=1, kl_scaling=None, *args, **kwargs):
    """Initialize inference algorithm. It initializes hyperparameters
    and builds ops for the algorithm's computation graph.
    Args:
      n_samples: int.
        Number of samples from variational model for calculating
        stochastic gradients.
      kl_scaling: dict of RandomVariable to tf.Tensor.
        Provides option to scale terms when using ELBO with KL divergence.
        If the KL divergence terms are
        $\\alpha_p \mathbb{E}_{q(z\mid x, \lambda)} [
              \log q(z\mid x, \lambda) - \log p(z)],$
        then pass {$p(z)$: $\\alpha_p$} as `kl_scaling`,
        where $\\alpha_p$ is a tensor. Its shape must be broadcastable;
        it is multiplied element-wise to the batchwise KL terms.
    """
    if kl_scaling is None:
      kl_scaling = {}
    if n_samples <= 0:
      raise ValueError(
          "n_samples should be greater than zero: {}".format(n_samples))

    self.n_samples = n_samples
    self.kl_scaling = kl_scaling
    return super(test, self).initialize(*args, **kwargs)

  def build_loss_and_gradients(self, var_list):
    """Wrapper for the `KLqp` loss function.
    $-\\text{ELBO} =
        -\mathbb{E}_{q(z; \lambda)} [ \log p(x, z) - \log q(z; \lambda) ]$
    KLqp supports
    1. score function gradients [@paisley2012variational]
    2. reparameterization gradients [@kingma2014auto]
    of the loss function.
    If the KL divergence between the variational model and the prior
    is tractable, then the loss function can be written as
    $-\mathbb{E}_{q(z; \lambda)}[\log p(x \mid z)] +
        \\text{KL}( q(z; \lambda) \| p(z) ),$
    where the KL term is computed analytically [@kingma2014auto]. We
    compute this automatically when $p(z)$ and $q(z; \lambda)$ are
    Normal.
    """
    is_reparameterizable = all([
        rv.reparameterization_type ==
        tf.contrib.distributions.FULLY_REPARAMETERIZED
        for rv in six.itervalues(self.latent_vars)])
    is_analytic_kl = all([isinstance(z, Normal) and isinstance(qz, Normal)
                          for z, qz in six.iteritems(self.latent_vars)])
    if not is_analytic_kl and self.kl_scaling:
      raise TypeError("kl_scaling must be None when using non-analytic KL term")
    if is_reparameterizable:
      if is_analytic_kl:
        return build_reparam_kl_loss_and_gradients(self, var_list)
      # elif is_analytic_entropy:
      #    return build_reparam_entropy_loss_and_gradients(self, var_list)
      else:
        return build_reparam_loss_and_gradients(self, var_list)
    else:
      # Prefer Rao-Blackwellization over analytic KL. Unknown what
      # would happen stability-wise if the two are combined.
      # if is_analytic_kl:
      #   return build_score_kl_loss_and_gradients(self, var_list)
      # Analytic entropies may lead to problems around
      # convergence; for now it is deactivated.
      # elif is_analytic_entropy:
      #    return build_score_entropy_loss_and_gradients(self, var_list)
      # else:
      return build_score_rb_loss_and_gradients(self, var_list)


def build_reparam_kl_loss_and_gradients(inference, var_list):
  """Build loss function. Its automatic differentiation
  is a stochastic gradient of
  .. math::
    -\\text{ELBO} =  - ( \mathbb{E}_{q(z; \lambda)} [ \log p(x \mid z) ]
          + \\text{KL}(q(z; \lambda) \| p(z)) )
  based on the reparameterization trick [@kingma2014auto].
  It assumes the KL is analytic.
  Computed by sampling from $q(z;\lambda)$ and evaluating the
  expectation using Monte Carlo sampling.
  """
 
  p_log_lik = [0.0] * inference.n_samples
  base_scope = tf.get_default_graph().unique_name("inference") + '/'
  for s in range(inference.n_samples):
    # Form dictionary in order to replace conditioning on prior or
    # observed variable with conditioning on a specific value.
    scope = base_scope + tf.get_default_graph().unique_name("sample")
    dict_swap = {}
    for x, qx in six.iteritems(inference.data):
      if isinstance(x, RandomVariable):
        if isinstance(qx, RandomVariable):
          qx_copy = copy(qx, scope=scope)
          dict_swap[x] = qx_copy.value()
        else:
          dict_swap[x] = qx

    for z, qz in six.iteritems(inference.latent_vars):
      # Copy q(z) to obtain new set of posterior samples.
      qz_copy = copy(qz, scope=scope)
      dict_swap[z] = qz_copy.value()

    for x in six.iterkeys(inference.data):
      if isinstance(x, RandomVariable):
        x_copy = copy(x, dict_swap, scope=scope)
        p_log_lik[s] += tf.reduce_sum(
            inference.scale.get(x, 1.0) * x_copy.log_prob(dict_swap[x]))

  p_log_lik = tf.reduce_mean(p_log_lik)

  kl_penalty = tf.reduce_sum([
      tf.reduce_sum(inference.kl_scaling.get(z, 1.0) * kl_divergence(qz, z))
      for z, qz in six.iteritems(inference.latent_vars)])

  if inference.logging:
    tf.summary.scalar("loss/p_log_lik", p_log_lik,
                      collections=[inference._summary_key])
    tf.summary.scalar("loss/kl_penalty", kl_penalty,
                      collections=[inference._summary_key])

  loss = -(p_log_lik - kl_penalty)
  grads = tf.gradients(loss, var_list)
  grads_and_vars = list(zip(grads, var_list))
  return loss, grads_and_vars


# def build_reparam_kl_loss_and_gradients(inference, var_list):
#   """Build loss function. Its automatic differentiation
#   is a stochastic gradient of
#   .. math::
#     -\\text{ELBO} =  - ( \mathbb{E}_{q(z; \lambda)} [ \log p(x \mid z) ]
#           + \\text{KL}(q(z; \lambda) \| p(z)) )
#   based on the reparameterization trick [@kingma2014auto].
#   It assumes the KL is analytic.
#   Computed by sampling from $q(z;\lambda)$ and evaluating the
#   expectation using Monte Carlo sampling.
#   """
 
#   p_log_lik = [0.0] * inference.n_samples
#   base_scope = tf.get_default_graph().unique_name("inference") + '/'
#   for s in range(inference.n_samples):
#     # Form dictionary in order to replace conditioning on prior or
#     # observed variable with conditioning on a specific value.
#     scope = base_scope + tf.get_default_graph().unique_name("sample")
#     dict_swap = {}
#     for x, qx in six.iteritems(inference.data):
#       if isinstance(x, RandomVariable):
#         if isinstance(qx, RandomVariable):
#           qx_copy = copy(qx, scope=scope)
#           dict_swap[x] = qx_copy.value()
#         else:
#           dict_swap[x] = qx

#     for z, qz in six.iteritems(inference.latent_vars):
#       # Copy q(z) to obtain new set of posterior samples.
#       qz_copy = copy(qz, scope=scope)
#       dict_swap[z] = qz_copy.value()

#     for x in six.iterkeys(inference.data):
#       if isinstance(x, RandomVariable):
#         x_copy = copy(x, dict_swap, scope=scope)
#         p_log_lik[s] += tf.reduce_sum(
#             inference.scale.get(x, 1.0) * x_copy.log_prob(dict_swap[x]))

#   p_log_lik = tf.reduce_mean(p_log_lik)

#   kl_penalty = tf.reduce_sum([
#       tf.reduce_sum(inference.kl_scaling.get(z, 1.0) * kl_divergence(qz, z))
#       for z, qz in six.iteritems(inference.latent_vars)])

#   reg_penalty = tf.reduce_sum(tf.losses.get_regularization_losses())

#   if inference.logging:
#     tf.summary.scalar("loss/p_log_lik", p_log_lik,
#                       collections=[inference._summary_key])
#     tf.summary.scalar("loss/kl_penalty", kl_penalty,
#                       collections=[inference._summary_key])
#     tf.summary.scalar("loss/reg_penalty", reg_penalty,
#                       collections=[inference._summary_key])

#   loss = -(p_log_lik - kl_penalty - reg_penalty)
#   grads = tf.gradients(loss, var_list)
#   grads_and_vars = list(zip(grads, var_list))
#   return loss, grads_and_vars
