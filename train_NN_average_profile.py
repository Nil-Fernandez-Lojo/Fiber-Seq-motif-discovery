import numpy as np
from util import load_genome, split_average_profile_with_peaks_to_fixed_size_inputs,load_bed_file_peaks,generator_data_average_profile,get_mu_and_std_average_profile_peaks,generator_input_average_profile,predictions_NN_to_wig_file
from NN_functions import Model, train_on_average_profile
import tensorflow as tf
import wandb
import os
os.environ["CUDA_VISIBLE_DEVICES"]="0" # second gpu

np.random.seed(0)
tf.random.set_seed(1)
genome_file = "dm6.fa"
average_profile_file_npy = "average_profile/profile_smoothed_w_2.npy"
bed_file_peaks = 'methyl_bed_file/GSM4411217_DS68346_hotspot_fdr05.bed'
train_predictions_file = 'train_predictions.wig'
val_predictions_file = 'validation_predictions.wig'

file_out_model = 'model_trained'

configs = {"learning_rate": 0.001,
			"epochs": 10,
			"batch_size": 128,
			"half_width_output": 500, #Output length is 2*half_width_output+1
			"architecture": "BPnet",
			"loss_fn": "RMSE",
			"n_layers": 5,
			"input_augmentation":1000,
			"dataset": average_profile_file_npy,
			"n_units_layer": 64,
			"filter_1_size": 25,
			"filter_other_size":3
			}

#### NN architecture

run = wandb.init(project="fiber-training", entity="lcb_2", config=configs)
config = wandb.config

w = config.half_width_output
input_augmentation = config.input_augmentation
n_units_layer = config.n_units_layer
filter_1_size = config.filter_1_size
filter_other_size = config.filter_other_size

ft = 0.7
fv = 0.2

fraction_test = 1-ft-fv

output_length = 2*w+1
input_length = output_length + 2*input_augmentation

chr_names,chr_seq = load_genome(genome_file)
profile = np.load(average_profile_file_npy)
peaks = load_bed_file_peaks(bed_file_peaks,chr_names,chr_seq)

x,y = split_average_profile_with_peaks_to_fixed_size_inputs(profile,chr_seq,peaks,w)
x_train_before_shuffling = x[:round(ft*len(x)), :]
y_train = y[:round(ft*len(x)), :]
x_val = x[round(ft*len(x)):round((ft+fv)*len(x)), :]
y_val = y[round(ft*len(x)):round((ft+fv)*len(x)), :]

p = np.random.permutation(len(x_train_before_shuffling))
x_train = x_train_before_shuffling[p,:]
y_train = y_train[p,:]


model = Model(input_length, output_length,n_units_layer,config.n_layers, filter_1_size,filter_other_size, average_profile = True)
if config.loss_fn =="KL":
	loss_fn = tf.keras.losses.KLDivergence()
elif config.loss_fn =="cross entropy":
	loss_fn = tf.keras.losses.BinaryCrossentropy()
elif config.loss_fn =="RMSE":
	loss_fn = tf.keras.losses.MeanSquaredError()
else:
	raise("Only avaliable loss functions are for now KL, cross entropy and RMSE")
optimizer = tf.keras.optimizers.Adam(learning_rate=config.learning_rate)

train_acc_metric = tf.keras.metrics.RootMeanSquaredError()
val_acc_metric = tf.keras.metrics.RootMeanSquaredError()

mu_y_train, std_y_train = get_mu_and_std_average_profile_peaks(y_train)
mu_y_val, std_y_val = get_mu_and_std_average_profile_peaks(y_val)
print('training profile mean: ', mu_y_train, ' std: ', std_y_train)
print('validation profile mean: ', mu_y_val, ' std: ', std_y_val)

#model.compile(optimizer =optimizer, loss = loss_fn )
#model.fit(generator_data_average_profile(x_train,y_train,chr_seq,input_augmentation,config.batch_size),
#	epochs = config.epochs,
#	verbose = 2,
#	steps_per_epoch = round(len(x_train)/config.batch_size))
train_on_average_profile(model,
		optimizer,
		loss_fn,
		config.epochs,
		config.batch_size,
		x_train,
		y_train,
		x_val,
		y_val,
		chr_seq,
		input_augmentation,
		train_acc_metric,
		val_acc_metric)
model.save(file_out_model)

y_train_pred = model.predict(generator_input_average_profile(x_train_before_shuffling,
															chr_seq,
															input_augmentation,
															config.batch_size),
							steps = np.floor(len(x_train_before_shuffling)/config.batch_size)
)

y_val_pred = model.predict(generator_input_average_profile(x_val,
															chr_seq,
															input_augmentation,
															config.batch_size),
							steps = np.floor(len(x_val)/config.batch_size)
)
print("shape x_train_before_shuffling", x_train_before_shuffling.shape)
print("shape y_train_pred", y_train_pred.shape)
print("shape x_val", x_val.shape)
print("shape y_val_pred", y_val_pred.shape)

predictions_NN_to_wig_file(x_train_before_shuffling,y_train_pred,chr_names, train_predictions_file)
predictions_NN_to_wig_file(x_val,y_val_pred,chr_names, val_predictions_file)
run.join()

