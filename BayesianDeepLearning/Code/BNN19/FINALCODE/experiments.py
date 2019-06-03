import tensorflow as tf
import numpy as np

import tensorflow_probability as tfp
tfd = tfp.distributions

import warnings
warnings.filterwarnings('ignore')

from tensorflow.examples.tutorials.mnist import input_data
mnist = input_data.read_data_sets('MNIST_DATA', one_hot=True)

from convnet import ConvNet


def run_experiment(algo,average_gradients, batch_size, iterations, verbose):
  batch_size = batch_size
  tf.reset_default_graph()
  net = ConvNet()

  validation_batch = mnist.test.images
  val_count = validation_batch.shape[0]
  validation_batch = np.reshape(validation_batch, (val_count, 28, 28, 1))
  validation_labels = mnist.test.labels

  net.setup_train(algo,average_gradients=average_gradients)
  training_log = []
  listloss= []
  with tf.Session() as sess:
    sess.run(tf.global_variables_initializer())
    for i in range(iterations):
      batch = mnist.train.next_batch(batch_size)
      input_batch = np.reshape(batch[0], (batch_size, 28, 28, 1))
      loss = net.train(sess, input_batch, batch[1])
      listloss.append(loss)
      if (i+1) % 100 == 0:
        accuracy = net.evaluate(sess, validation_batch, validation_labels)
        training_log.append((accuracy, i+1))
        if verbose:
          print('[{:d}/{:d}] loss: {:.3g}, accuracy: {:.3g}%'.format(i+1, iterations, loss, accuracy))
    accuracy = net.evaluate(sess, validation_batch, validation_labels)
    training_log.append((accuracy, iterations))
    best = sorted(training_log, key=lambda x: x[0], reverse=True)[0]
    print('Training finished. Best accuracy: {:.3g} at iteration {:d}.'.format(best[0], best[1]))
    return listloss


batch_size = 128
iterations = 100

adam = run_experiment(algo='adam', average_gradients=1, batch_size = batch_size, iterations=iterations, verbose= True)
bbb = run_experiment(algo='bbb', average_gradients=1, batch_size = batch_size, iterations=iterations, verbose= True)
momentum = run_experiment(algo='momentum', average_gradients=1, batch_size = batch_size, iterations=iterations, verbose= True)


###PLOTS
import matplotlib.pyplot as plt
iter = len(adam)
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16, 3.5))
plt.plot(np.arange(iter), adam , label='ADAM')
plt.plot(np.arange(iter), bbb , label='BBB')
#plt.plot(np.arange(iter), momentum , label='Momentum')
leg = plt.legend(fontsize=20,fancybox=True, loc='right')
leg.get_frame().set_alpha(0.5)
plt.xlabel('Epoch', fontsize=15)
plt.ylabel('Negated ELBO', fontsize=15)
plt.show()
