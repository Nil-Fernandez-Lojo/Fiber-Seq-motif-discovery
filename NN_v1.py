#import tensorflow as tf
from util import load_genome,load_bed_file,break_reads_to_fixed_size,shuffle_data
import time

#### Parameters NN
output_length = 1000
input_augmentation = 1000
n_units_layer = 64
filter_1_size = 25
filter_other_size = 3

input_length = output_length + 2*input_augmentation


chr_names,chr_seq = load_genome("dm6.fa")
seq_pos,n_methyl,methyl_rel_pos = load_bed_file("GSM4411218_tracks_m6A_DS75167_10.dm6.bed",chr_names)
seq_pos,n_methyl,methyl_rel_pos = break_reads_to_fixed_size(seq_pos,n_methyl,methyl_rel_pos,output_length)

#check umbalanced training set and loss function
#split training and validation


# #Still need to add direct connexion input and output and fix weight
# model = tf.keras.Sequential()
# model.add(tf.keras.layers.Conv1D(n_units_layer, kernel_size=filter_1_size,padding='same', activation='relu', input_shape=(input_length, 4)))

# for i in range(1, 10):
# 	model.add(tf.keras.layers.Conv1D(n_units_layer, kernel_size=filter_other_size, padding='same',activation='relu', dilation_rate=2**i))
# model.add(tf.keras.layers.Dense(output_length,activation="sigmoid"))

# loss_fn = tf.keras.losses.BinaryCrossentropy()
# model.compile(optimizer='adam',
#               loss=loss_fn)

# model.fit_generator(x_train, y_train, epochs=5) #check of this works and train with like 5% of the samples
#model.evaluate(x_test,  y_test, verbose=2)
