
# coding: utf-8

# In[167]:

from __future__ import print_function
import numpy as np
import tensorflow as tf
from six.moves import cPickle as pickle
from six.moves import range
import math


# In[168]:

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


        
dataset = []
for d in data:
    d1 = []
    for l in d:
        d1.append(np.float32(ord(l)))
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



# In[215]:



shuf = np.append(dataset,labels,axis = 1)
np.random.shuffle(shuf)

dataset = shuf[:,:20]
labels = shuf[:,20:22]
print (np.shape(dataset),np.shape(labels))


# In[216]:


    
trainnum = 200000
testnum = 20000
train_dataset = dataset[:trainnum]
train_labels = labels[:trainnum]

test_dataset = dataset[trainnum:trainnum + testnum]
test_labels = labels[trainnum:trainnum + testnum]

valid_dataset = dataset[trainnum + testnum:trainnum + 2*testnum]
valid_labels = labels[trainnum + testnum:trainnum + 2*testnum]



# In[217]:

print (dataset[5],)


# In[218]:



graph = tf.Graph()
with graph.as_default():

  # Input data.
  
  tf_train_dataset = tf.constant(train_dataset)
  tf_train_labels = tf.constant(train_labels)
  tf_valid_dataset = tf.constant(valid_dataset)
  tf_test_dataset = tf.constant(test_dataset)
  
  # Variables.
 
  weights = tf.Variable(
    tf.truncated_normal([20, 2]))
  biases = tf.Variable(tf.zeros([2]))
  
  # Training computation.
  
  logits = tf.matmul(tf_train_dataset, weights) + biases
  loss = tf.reduce_mean(
    tf.nn.softmax_cross_entropy_with_logits(logits, tf_train_labels))
  
  # Optimizer.

  optimizer = tf.train.GradientDescentOptimizer(0.01).minimize(loss)
  
  # Predictions
  train_prediction = tf.nn.softmax(logits)
  valid_prediction = tf.nn.softmax(
    tf.matmul(tf_valid_dataset, weights) + biases)
  test_prediction = tf.nn.softmax(tf.matmul(tf_test_dataset, weights) + biases)


# In[ ]:




# In[219]:

num_steps = 1601

def accuracy(predictions, labels):
  return (100.0 * np.sum(np.argmax(predictions, 1) == np.argmax(labels, 1))
          / predictions.shape[0])

with tf.Session(graph=graph) as session:
   
  tf.initialize_all_variables().run()
  print('Initialized')
  for step in range(num_steps):
    
    _, l, predictions = session.run([optimizer, loss, train_prediction])
    if (step % 100 == 0):
      print('Loss at step %d: %f' % (step, l))
      print('Training accuracy: %.1f%%' % accuracy(
        predictions, train_labels))
      
      print('Validation accuracy: %.1f%%' % accuracy(
        valid_prediction.eval(), valid_labels))
  print('Test accuracy: %.1f%%' % accuracy(test_prediction.eval(), test_labels))


# In[ ]:



