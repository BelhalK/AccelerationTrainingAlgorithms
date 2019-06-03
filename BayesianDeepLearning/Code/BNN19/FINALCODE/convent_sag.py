import tensorflow as tf
import numpy as np

import tensorflow_probability as tfp
tfd = tfp.distributions

import warnings
warnings.filterwarnings('ignore')

from tensorflow.examples.tutorials.mnist import input_data
mnist = input_data.read_data_sets('MNIST_DATA', one_hot=True)

from convnet import ConvNet


class ConvNet(object):
  def __init__(self, dropout=0):
    #self._input = tf.placeholder(tf.float32, shape=[None, 28, 28, 1])
    self._input = tf.keras.layers.Input(shape=[28, 28, 1], dtype='float32')
    self._training = tf.placeholder_with_default(False, shape=[])
    self._labels = tf.placeholder(tf.float32, shape=[None, 10])
    

    if algo=="dropout":
      dropout=0.5

    out = self._input 
    if dropout:
      out = tf.layers.dropout(out, rate=0.2 if dropout > 0.2 else dropout, training=self._training)
    out = tfp.layers.Convolution2DReparameterization(6,kernel_size=5,padding="SAME",activation=tf.nn.relu)(out)
    out = tf.layers.max_pooling2d(out, pool_size=2, strides=2) # 12
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
    
    optimizer = tf.train.GradientDescentOptimizer(learning_rate=lr)

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
 


def run_experiment(average_gradients, batch_size, iterations, verbose):
  batch_size = batch_size
  tf.reset_default_graph()
  net = ConvNetMISSO()

  net.setup_train(average_gradients=average_gradients)
  training_log = []
  listloss= []
  with tf.Session() as sess:
    sess.run(tf.global_variables_initializer())
    indivgrads = []
    for indiv in range(0,total,batch_size):
      print(indiv)
      grads = tf.gradients(net._loss_op, tf.trainable_variables())
      var_updates = []
      var_list = tf.trainable_variables()
      for grad, var in zip(grads, var_list):
          var_updates.append(var.assign_sub(0.001 * grad))
      net._train_op = tf.group(*var_updates)
      indivgrads.append(grads)
    for i in range(iterations):
      print(i)
      batch = mnist.train.next_batch(batch_size)
      input_batch = np.reshape(batch[0], (batch_size, 28, 28, 1))
      loss = net.train(sess, input_batch, batch[1])
      listloss.append(loss)
      if (i+1) % 100 == 0:
        print('[{:d}/{:d}] loss: {:.3g}, accuracy: {:.3g}%'.format(i+1, iterations, loss, accuracy))
    print('Training finished.')
    return listloss



def run_experimentsag(average_gradients, batch_size, iterations, verbose):
  batch_size = batch_size
  tf.reset_default_graph()
  net = ConvNetMISSO()

  net.setup_train(average_gradients=average_gradients)
  training_log = []
  listloss= []
  with tf.Session() as sess:
    sess.run(tf.global_variables_initializer())
    indivgrads = []
    for indiv in range(0,total,batch_size):
      print(indiv)
      grads = tf.gradients(net._loss_op, tf.trainable_variables())
      var_updates = []
      var_list = tf.trainable_variables()
      for grad, var in zip(grads, var_list):
          var_updates.append(var.assign_sub(0.001 * grad))
      net._train_op = tf.group(*var_updates)
      indivgrads.append(grads)
    for epoch in range(iterations):
      for index in range(0,int(total/batch_size)):
        print(index)
        batch = mnist.train.next_batch(batch_size)
        input_batch = np.reshape(batch[0], (batch_size, 28, 28, 1))
        feed_dict = {net._input: input_batch,net._labels: batch[1],net._training: True}

        grads = tf.gradients(net._loss_op, tf.trainable_variables())
        indivgrads[index] = grads
        var_updates = []
        var_list = tf.trainable_variables()
        print('ok')
        for gradstemp in indivgrads:
            for grad, var in zip(gradstemp, var_list):
                var_updates.append(var.assign_sub(0.001 * grad)) #\theta^k - \grad f_{\tau_i^k}
        net._train_op = tf.group(*var_updates)

        sess.run(net._train_op,feed_dict=feed_dict)
        loss = sess.run(net._loss_op, feed_dict=feed_dict)
        listloss.append(loss)
    print('Training finished.')
    return listloss


batch_size = 128
iterations = 100

sag = run_experimentsag(average_gradients=1, batch_size = batch_size, iterations=iterations, verbose= True)
