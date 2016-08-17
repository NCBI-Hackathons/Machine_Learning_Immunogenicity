
# coding: utf-8

# In[1]:

from __future__ import print_function
import numpy as np
import tensorflow as tf
from six.moves import cPickle as pickle
from six.moves import range
import math


# In[3]:

#load data,channge type
load = np.genfromtxt("/Users/yueyuchen/Documents/ML/ML_Immune/Machine_Learning_Immunogenicity/data/bcell.csv",dtype = str, delimiter = ',')

length1 = np.shape(load)[0]
labels1 = np.zeros(length1)
load1 = np.zeros(length1)
for i in range(length1):
    labels1[i] = (load[i,1] == 'Positive')
for i in range(length1):   
    load1[i] = len(load[i,0])
data = []
labels = []
for i in range(length1):
    if(load1[i] < 20):
        data.append(load[i,0])
        labels.append(labels1[i])

length = len(data)

Alph = {'A':1,'a':1,
        'B':2,'b':2,
        'C':3,'c':3,
        'D':4,'d':4,
        'E':5,'e':5,
        'F':6,'f':6,
        'G':7,'g':7,
        'H':8,'h':8,
        'I':9,'i':9,
        'J':10,'j':10,
        'K':11,'k':11,
        'L':12,'l':12,
        'M':13,'m':13,
        'N':14,'n':14,
        'O':15,'o':15,
        'P':16,'p':16,
        'Q':17,'q':17,
        'R':18,'r':18,
        'S':19,'s':19,
        'T':20,'t':20,
        'U':21,'u':21,
        'V':22,'v':22,
        'W':23,'w':23,
        'X':24,'x':24,
        'Y':25,'y':25,
        'Z':26,'z':26
        }
        
dataset = []
for d in data:
    d1 = []
    for l in d:
        d1.append(np.float32(Alph[l]))
    for j in range(20-len(d1)):
        d1.append(np.float32(0))
    dataset.append(d1)

for i in range(len(labels)):
    if (labels[i] == 0):
        labels[i] = [1,0]
    else:
        labels[i] = [0,1]
        
for i in range(length):
    dataset[i] = np.asarray(dataset[i])
    labels[i] = np.float32(np.asarray(labels[i]))

dataset = np.asarray(dataset)
labels = np.asarray(labels)


# In[12]:

#shuffle dataset
shuf = np.append(dataset,labels,axis = 1)
np.random.shuffle(shuf)

dataset = shuf[:,:20]
labels = shuf[:,20:22]
print (np.shape(dataset),np.shape(labels))


# In[13]:

trainnum = 200000
testnum = 20000
train_dataset = dataset[:trainnum]
train_labels = labels[:trainnum]

test_dataset = dataset[trainnum:trainnum + testnum]
test_labels = labels[trainnum:trainnum + testnum]

valid_dataset = dataset[trainnum + testnum:trainnum + 2*testnum]
valid_labels = labels[trainnum + testnum:trainnum + 2*testnum]


# In[14]:

batch_size = 1000
num_hidden_nodes = 500
peplength = 20
num_labels = 2

graph = tf.Graph()
with graph.as_default():

  # Input data.
  tf_train_dataset = tf.placeholder(tf.float32,
                                    shape=(batch_size, peplength))
  tf_train_labels = tf.placeholder(tf.float32, shape=(batch_size, num_labels))
  tf_valid_dataset = tf.constant(valid_dataset)
  tf_test_dataset = tf.constant(test_dataset)
  
  # Variables.
  weights1 = tf.Variable(
    tf.truncated_normal([peplength, num_hidden_nodes]))
  biases1 = tf.Variable(tf.zeros([num_hidden_nodes]))
  weights2 = tf.Variable(
    tf.truncated_normal([num_hidden_nodes, num_labels]))
  biases2 = tf.Variable(tf.zeros([num_labels]))
  
  # Training computation.
  lay1_train = tf.nn.relu(tf.matmul(tf_train_dataset, weights1) + biases1)
  logits = tf.matmul(lay1_train, weights2) + biases2
  loss = tf.reduce_mean(
    tf.nn.softmax_cross_entropy_with_logits(logits, tf_train_labels))
  
  # Optimizer.
  optimizer = tf.train.GradientDescentOptimizer(0.0001).minimize(loss)
  
  # Predictions for the training, validation, and test data.
  train_prediction = tf.nn.softmax(logits)
  lay1_valid = tf.nn.relu(tf.matmul(tf_valid_dataset, weights1) + biases1)
  valid_prediction = tf.nn.softmax(tf.matmul(lay1_valid, weights2) + biases2)
  lay1_test = tf.nn.relu(tf.matmul(tf_test_dataset, weights1) + biases1)
  test_prediction = tf.nn.softmax(tf.matmul(lay1_test, weights2) + biases2)


# In[15]:

def accuracy(predictions, labels):
  return (100.0 * np.sum(np.argmax(predictions, 1) == np.argmax(labels, 1))
          / predictions.shape[0])


# In[16]:

num_steps = 15001

with tf.Session(graph=graph) as session:
  tf.initialize_all_variables().run()
  print("Initialized")
  for step in range(num_steps):
    
    offset = (step * batch_size) % (train_labels.shape[0] - batch_size)
    # Generate a minibatch.
    batch_data = train_dataset[offset:(offset + batch_size), :]
    batch_labels = train_labels[offset:(offset + batch_size), :]
    
    feed_dict = {tf_train_dataset : batch_data, tf_train_labels : batch_labels}
    _, l, predictions = session.run(
      [optimizer, loss, train_prediction], feed_dict=feed_dict)
    if (step % 500 == 0):
      print("Minibatch loss at step %d: %f" % (step, l))
      print("Minibatch accuracy: %.1f%%" % accuracy(predictions, batch_labels))
      print("Validation accuracy: %.1f%%" % accuracy(
        valid_prediction.eval(), valid_labels))
  print("Test accuracy: %.1f%%" % accuracy(test_prediction.eval(), test_labels))


# In[ ]:



