import numpy as np
import tensorflow as tf
from util import load_genome,load_bed_file,break_reads_to_fixed_size,shuffle_data,split_training_testing,generator_data
import time
import wandb
from wandb.keras import WandbCallback

np.random.seed(0)
wandb.init(project="fiber-training", entity="lcb_2")

#### Parameters NN
genome_file = "dm6.fa"
bed_file = "GSM4411218_tracks_m6A_DS75167_100.dm6.bed"
output_length = 1000
input_augmentation = 1000
n_units_layer = 64
filter_1_size = 25
filter_other_size = 3

batch_size = 32

fraction_training = 0.7

input_length = output_length + 2*input_augmentation

start = time.time()
chr_names,chr_seq = load_genome(genome_file)
seq_pos,n_methyl,methyl_rel_pos = load_bed_file(bed_file,chr_names)
seq_pos,n_methyl,methyl_rel_pos = break_reads_to_fixed_size(seq_pos,n_methyl,methyl_rel_pos,output_length)
seq_pos,n_methyl,methyl_rel_pos = shuffle_data(seq_pos,n_methyl,methyl_rel_pos)
seq_pos_training,n_methyl_training,methyl_rel_pos_training, seq_pos_testing,n_methyl_testing,methyl_rel_pos_testing = split_training_testing(seq_pos,n_methyl,methyl_rel_pos,fraction_training)
n_training_samples = len(seq_pos_training)
print(time.time()- start)

#Still need to add direct connexion input and output and fix weight
model = tf.keras.Sequential()
model.add(tf.keras.layers.Conv1D(n_units_layer, kernel_size=filter_1_size,padding='same', activation='relu', input_shape=(input_length, 4)))

#for i in range(1, 10):
#	model.add(tf.keras.layers.Conv1D(n_units_layer, kernel_size=filter_other_size, padding='same',activation='relu', dilation_rate=2**i))
model.add(tf.keras.layers.Flatten())
model.add(tf.keras.layers.Dense(output_length,activation="sigmoid"))

loss_fn = tf.keras.losses.BinaryCrossentropy()
model.compile(optimizer='adam',
              loss=loss_fn)


model.fit(x=generator_data(seq_pos_training,n_methyl_training,methyl_rel_pos_training, chr_seq,input_augmentation,batch_size),
	steps_per_epoch = round(n_training_samples/batch_size),
	epochs=5,
	verbose = 2) 
model.evaluate(generator_data(seq_pos_testing,n_methyl_testing,methyl_rel_pos_testing, chr_seq,input_augmentation,len(seq_pos_testing)) ,verbose=2)
