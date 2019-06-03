import tensorflow as tf
import numpy as np
import tensorflow_probability as tfp
tfd = tfp.distributions


class ConvNet(object):
  def __init__(self, dropout=0):
    #self._input = tf.placeholder(tf.float32, shape=[None, 28, 28, 1])
    self._input = tf.keras.layers.Input(shape=[28, 28, 1], dtype='float32')
    self._training = tf.placeholder_with_default(False, shape=[])
    self._labels = tf.placeholder(tf.float32, shape=[None, 10])
    
    
    out = self._input 
    if dropout:
      out = tf.layers.dropout(out, rate=0.2 if dropout > 0.2 else dropout, training=self._training)
    out = tfp.layers.Convolution2DReparameterization(6,kernel_size=5,padding="SAME",activation=tf.nn.relu)(out)
    if dropout:
      out = tf.layers.dropout(out, rate=dropout, training=self._training)
    out = tfp.layers.Convolution2DReparameterization(16,kernel_size=5,padding="SAME",activation=tf.nn.relu, name='layer1')(out)
    out = tf.layers.max_pooling2d(out, pool_size=2, strides=2) # 4
    if dropout:
      out = tf.layers.dropout(out, rate=dropout, training=self._training)
    out = tf.contrib.layers.flatten(out)
    out = tfp.layers.DenseReparameterization(120, activation=tf.nn.relu, name='layer2')(out)
    if dropout:
      out = tf.layers.dropout(out, rate=dropout, training=self._training)
    out = tfp.layers.DenseReparameterization(84, activation=tf.nn.relu, name='layer3')(out)
    if dropout:
      out = tf.layers.dropout(out, rate=dropout, training=self._training)
    out = tfp.layers.DenseReparameterization(10, name='layer4')(out)
    
    
    self._inference_op = out

    correct = tf.equal(tf.argmax(self._labels, 1), tf.argmax(self._inference_op, 1))
    self._accuracy_op = tf.reduce_mean(tf.cast(correct, tf.float32))

  @classmethod
  def conv2d(cls, input, filters, name=None):
    return tf.layers.conv2d(input, filters=filters, kernel_size=5, activation=tf.nn.relu, name=name)
  @classmethod
  def dense(cls, input, filters, name=None):
    return tf.layers.dense(input, filters, activation=tf.nn.relu, name=name)

  def setup_train(self,algo, average_gradients=1, lr=1e-3):
    self._average_gradients = average_gradients
    
    ###BNN loss
    labels_distribution = tfd.Categorical(logits=self._inference_op)
    # neg_log_likelihood = -tf.reduce_mean(input_tensor=labels_distribution.log_prob(self._labels))
    neg_log_likelihood = tf.losses.softmax_cross_entropy(onehot_labels=self._labels, logits=self._inference_op)
    model = tf.keras.Model(inputs=self._input, outputs=self._inference_op)
    
    kl = sum(model.losses) / 55000
    elbo_loss = neg_log_likelihood + kl
    self._loss_op = elbo_loss
    
    #CNN loss
    #self._loss_op = tf.losses.softmax_cross_entropy(onehot_labels=self._labels, logits=self._inference_op)
    
    if algo=="adam":
      optimizer = tf.train.AdamOptimizer(learning_rate=lr)
    if algo=="bbb":
      optimizer = tf.train.GradientDescentOptimizer(learning_rate=lr)
    if algo=="misso":
      optimizer = tf.train.GradientDescentOptimizer(learning_rate=lr)
    if algo=="momentum":
      optimizer = tf.train.MomentumOptimizer(learning_rate=lr, momentum=5e-5)

    if average_gradients == 1:
      # This 'train_op' computes gradients and applies them in one step.
      self._train_op = optimizer.minimize(self._loss_op)
    else:
      # here 'train_op' only applies gradients passed via placeholders stored
      # in 'grads_placeholders. The gradient computation is done with 'grad_op'.
      grads_and_vars = optimizer.compute_gradients(self._loss_op)
      avg_grads_and_vars = []
      self._grad_placeholders = []
      for grad, var in grads_and_vars:
        grad_ph = tf.placeholder(grad.dtype, grad.shape)
        self._grad_placeholders.append(grad_ph)
        avg_grads_and_vars.append((grad_ph, var))
      self._grad_op = [x[0] for x in grads_and_vars]
      self._train_op = optimizer.apply_gradients(avg_grads_and_vars)
      self._gradients = [] # list to store gradients

  def train(self, session, input_batch, output_batch):
    feed_dict = {
      self._input: input_batch,
      self._labels: output_batch,
      self._training: True
    }
    if self._average_gradients == 1:
      loss, _ = session.run([self._loss_op, self._train_op], feed_dict=feed_dict)
    else:
      loss, grads = session.run([self._loss_op, self._grad_op], feed_dict=feed_dict)
      self._gradients.append(grads)
      if len(self._gradients) == self._average_gradients:
        for i, placeholder in enumerate(self._grad_placeholders):
          feed_dict[placeholder] = np.stack([g[i] for g in self._gradients], axis=0).mean(axis=0)
        session.run(self._train_op, feed_dict=feed_dict)
        self._gradients = []
    return loss

  def evaluate(self, session, input_batch, output_labels):
    feed_dict = {
      self._input: input_batch,
      self._labels: output_labels
    }
    return session.run(self._accuracy_op, feed_dict=feed_dict) * 100
