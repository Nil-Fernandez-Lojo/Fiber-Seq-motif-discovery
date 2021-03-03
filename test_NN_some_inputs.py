import numpy as np
from util import load_genome, split_average_profile_with_peaks_to_fixed_size_inputs,load_bed_file_peaks,generator_data_average_profile
import tensorflow as tf
import os
os.environ["CUDA_VISIBLE_DEVICES"]="-1" # second gpu

np.random.seed(0)
genome_file = "dm6.fa"
average_profile_file_npy = "average_profile/profile_smoothed_w_5.npy"
bed_file_peaks = 'methyl_bed_file/GSM4411217_DS68346_hotspot_fdr05.bed'
batch_size = 5
w = 50
input_augmentation = 0

ft = 0.7
fv = 0.2

output_length = 2*w+1
input_length = output_length + 2*input_augmentation

chr_names,chr_seq = load_genome(genome_file)
profile = np.load(average_profile_file_npy)
peaks = load_bed_file_peaks(bed_file_peaks,chr_names,chr_seq)

x,y = split_average_profile_with_peaks_to_fixed_size_inputs(profile,chr_seq,peaks,w)
x_train = x[:round(ft*len(x)), :]
y_train = y[:round(ft*len(x)), :]
x_val = x[round(ft*len(x)):round((ft+fv)*len(x)), :]
y_val = y[round(ft*len(x)):round((ft+fv)*len(x)), :]

for (x_t, y_t) in generator_data_average_profile(x_val,
																			y_val,
																			chr_seq,
																			input_augmentation,
																			batch_size):
	break



model = tf.keras.models.load_model('model_trained',compile=False)
y_pred = model.predict(x_t)
print(model.layers[0].get_weights())
print(model.layers[1].get_weights())
print(model.layers[2].get_weights())
for i in range(len(y_pred)):
	print('test : ', i)
	print('sequence: ', chr_seq[x_val[i,0]][x_val[i,1]:x_val[i,2]])
	print("true profile: \n",y_t[i,:])
	print("predicted profile: \n",y_pred[i,:])
print( model.predict(x_t*0))
