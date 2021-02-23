import numpy as np
import tensorflow as tf
from util import load_genome,load_bed_file,break_reads_to_fixed_size,shuffle_data,split_training_testing
from NN_functions import Model, train
import time
import wandb

np.random.seed(0)


genome_file = "dm6.fa"
bed_file = "GSM4411218_tracks_m6A_DS75167_100.dm6.bed"
configs = {"learning_rate": 0.01,
			"epochs": 5,
			"batch_size": 32,
			"architecture": "BPnet",
			"dataset": bed_file}

#### NN architecture
output_length = 500
input_augmentation = 250
n_layers = 10
n_units_layer = 64
filter_1_size = 25
filter_other_size = 3


fraction_training = 0.7


#run = wandb.init(project="fiber-training", entity="lcb_2", config=configs)
#config = wandb.config

input_length = output_length + 2*input_augmentation

chr_names,chr_seq = load_genome(genome_file)
seq_pos,n_methyl,methyl_rel_pos = load_bed_file(bed_file,chr_names,chr_seq)
seq_pos,n_methyl,methyl_rel_pos = break_reads_to_fixed_size(seq_pos,n_methyl,methyl_rel_pos,output_length)
seq_pos,n_methyl,methyl_rel_pos = shuffle_data(seq_pos,n_methyl,methyl_rel_pos)
seq_pos_training,n_methyl_training,methyl_rel_pos_training, seq_pos_val,n_methyl_val,methyl_rel_pos_val = split_training_testing(seq_pos,n_methyl,methyl_rel_pos,fraction_training)

print("Fraction of bases methylated in training data: ",n_methyl_training[-1]/(np.sum(seq_pos_training[:,2]-seq_pos_training[:,1])))

model = Model(input_length, output_length,n_units_layer, n_layers, filter_1_size,filter_other_size)
model.summary()
exit()
loss_fn = tf.keras.losses.BinaryCrossentropy()
optimizer = tf.keras.optimizers.Adam(learning_rate=config.learning_rate)

train_acc_metric = tf.keras.metrics.BinaryAccuracy()
val_acc_metric = tf.keras.metrics.BinaryAccuracy()

train(model,
		optimizer,
		loss_fn,
		config.epochs,
		config.batch_size,
		seq_pos_training,
		n_methyl_training,
		methyl_rel_pos_training,
		seq_pos_val,
		n_methyl_val,
		methyl_rel_pos_val,
		chr_seq,
		input_augmentation,
		train_acc_metric,
		val_acc_metric)

run.join()


#model.evaluate(generator_data(seq_pos_testing,n_methyl_testing,methyl_rel_pos_testing, chr_seq,input_augmentation,len(seq_pos_testing)) ,verbose=2)
