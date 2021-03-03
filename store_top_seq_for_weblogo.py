import random
import numpy as np
import tensorflow as tf
from util import seq_to_1_hot
folder_seq_weblgo = 'seq_filter_1_weblogo'
number_seq = 10**6
number_top_seq = 100 #used to generate weblogo
alphabet = "acgt"

def randseq(length):
	return ''.join(random.choice(alphabet) for i in range(length))
def write_fasta(list_seq, file_name):
	f=open(file_name, 'w')
	for i in range(len(list_seq)):
		f.write('>'+str(i)+'\n')
		f.write(list_seq[i]+'\n')
	f.close()


model = tf.keras.models.load_model('model_trained',compile=False)
filter_values = model.layers[1].get_weights()[0]
filter_bias_values = model.layers[1].get_weights()[1]

n_filter = len(filter_bias_values)
size_filter = filter_values.shape[0]

scores_seq = np.zeros((number_seq, n_filter))
list_seq = [None] * number_seq

for i in range(number_seq):
	if (i%1000 ==0): print(i,"out of",number_seq)
	list_seq[i] = randseq(size_filter)
	for j in range(n_filter):
		seq_1_hot = seq_to_1_hot(list_seq[i],0) #25x4
		scores_seq[i,j] = filter_bias_values[j] + np.sum(seq_1_hot * filter_values[:,:,j])

list_seq = np.array(list_seq)
scores_seq[scores_seq<0] = 0

top_seq_idx = np.zeros((number_top_seq, n_filter), dtype = int)
for j in range(n_filter):
	top_seq_idx[:,j] = np.argpartition(scores_seq[:,j], -number_top_seq)[-number_top_seq:].astype(int)
	top_seq_filter = list_seq[top_seq_idx[:,j]]
	write_fasta(top_seq_filter, folder_seq_weblgo+"/filter_"+str(j)+".fa")




