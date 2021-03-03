# Most of the functions come from https://colab.research.google.com/drive/126c1k5IfbQpE7dVmhnoDTdmFfC7CgJqg?usp=sharing#scrollTo=rdU2jDOB4itB
import wandb
from wandb.keras import WandbCallback
import tensorflow as tf
from util import generator_data,generator_data_average_profile
import numpy as np

def Model(input_length, output_length,n_units_layer, n_layers, filter_1_size,filter_other_size,average_profile = False,T_can_be_methyl = False):
	#Still need to add direct connexion input and output and fix weight
	#Careful, filter_1_size must be an odd number!
	s = int((input_length - output_length)/2)
	e = s + output_length

	inputs = tf.keras.Input(shape=(input_length, 4))
	x = tf.keras.layers.Conv1D(n_units_layer, 
		kernel_size=filter_1_size,
		padding='same', 
		activation='relu')(inputs)
	for i in range(1, n_layers):
		conv_x = tf.keras.layers.Conv1D(n_units_layer, 
			kernel_size=filter_other_size, 
			padding='same',
			activation='relu', 
			dilation_rate=2**i)(x)
		x = tf.keras.layers.add([conv_x, x])
	#x = tf.keras.layers.Flatten()(x)
	#x = tf.keras.layers.GlobalAveragePooling1D(data_format='channels_first')(x)
	#x = tf.keras.layers.Dense(output_length,activation="sigmoid")(x)
	
	x = tf.keras.layers.Reshape((-1, 1, n_units_layer))(x)
	x = tf.keras.layers.Conv2DTranspose(1, 
		kernel_size=(filter_1_size, 1), 
		padding='same')(x)
	x = tf.keras.layers.Flatten()(x)
	x = x[:,s:e]
	x = tf.keras.activations.sigmoid(x)

	if not average_profile:
		inputs_potentially_methyl = inputs[:,s:e,0]
		if T_can_be_methyl: inputs_potentially_methyl += inputs[:,s:e,1]
		x = tf.keras.layers.Multiply()([inputs_potentially_methyl, x])

	model = tf.keras.Model(inputs, x)

	# model = tf.keras.Sequential()
	# model.add(tf.keras.layers.Conv1D(n_units_layer, kernel_size=filter_1_size,padding='same', activation='relu', input_shape=(input_length, 4)))
	# for i in range(1, n_layers):
	# 	model.add(tf.keras.layers.Conv1D(n_units_layer, kernel_size=filter_other_size, padding='same',activation='relu', dilation_rate=2**i))

	# model.add(tf.keras.layers.Flatten())
	# model.add(tf.keras.layers.Dense(output_length,activation="sigmoid"))

	return model

def train_step(x, y, model, optimizer, loss_fn, train_acc_metric):
	with tf.GradientTape() as tape:
		logits = model(x, training=True)
		loss_value = loss_fn(y, logits)

	grads = tape.gradient(loss_value, model.trainable_weights)
	optimizer.apply_gradients(zip(grads, model.trainable_weights))

	train_acc_metric.update_state(y, logits)

	return loss_value

def test_step(x, y, model, loss_fn, val_acc_metric):
	val_logits = model(x, training=False)
	loss_value = loss_fn(y, val_logits)
	val_acc_metric.update_state(y, val_logits)

	return loss_value

#Need to change train_dataset val_dataset!!!

def train(model,
		optimizer,
		loss_fn,
		epochs,
		batch_size,
		seq_pos_training,
		n_methyl_training,
		methyl_rel_pos_training,
		seq_pos_val,
		n_methyl_val,
		methyl_rel_pos_val,
		chr_seq,
		input_augmentation,
		train_acc_metric,
		val_acc_metric):

	n_training_samples = len(seq_pos_training)
	training_step_per_epoch = round(n_training_samples/batch_size)

	n_val_samples = len(seq_pos_val)
	val_step_per_epoch = round(n_val_samples/batch_size)

	for epoch in range(epochs):
		print("\nStart of epoch %d" % (epoch,))

		train_loss = []   
		val_loss = []

		training_steps = 0
		# Iterate over the batches of the dataset
		for (x_batch_train, y_batch_train) in generator_data(seq_pos_training,
																	n_methyl_training,
																	methyl_rel_pos_training, 
																	chr_seq,
																	input_augmentation,
																	batch_size):
			loss_value = train_step(x_batch_train, y_batch_train, 
									model, optimizer, 
									loss_fn, train_acc_metric)
			train_loss.append(float(loss_value))
			training_steps+=1
			if (training_steps%training_step_per_epoch == 0): break


		# Run a validation loop at the end of each epoch
		val_steps = 0
		for (x_batch_val, y_batch_val) in generator_data(seq_pos_val,
																	n_methyl_val,
																	methyl_rel_pos_val, 
																	chr_seq,
																	input_augmentation,
																	batch_size):
			val_loss_value = test_step(x_batch_val, y_batch_val, 
									   model, loss_fn, 
									   val_acc_metric)
			val_loss.append(float(val_loss_value))

			val_steps+=1
			if (val_steps%val_step_per_epoch == 0): break
			
		# Display metrics at the end of each epoch
		train_acc = train_acc_metric.result()
		print("Training acc over epoch: %.4f" % (float(train_acc),))

		val_acc = val_acc_metric.result()
		print("Validation acc: %.4f" % (float(val_acc),))

		# Reset metrics at the end of each epoch
		train_acc_metric.reset_states()
		val_acc_metric.reset_states()

		# log metrics using wandb.log
		wandb.log({'epochs': epoch,
				   'loss': np.mean(train_loss),
				   'acc': float(train_acc), 
				   'val_loss': np.mean(val_loss),
				   'val_acc':float(val_acc)})

def train_on_average_profile(model,
		optimizer,
		loss_fn,
		epochs,
		batch_size,
		x_train,
		y_train,
		x_val,
		y_val,
		chr_seq,
		input_augmentation,
		train_acc_metric,
		val_acc_metric):

	n_training_samples = len(x_train)
	training_step_per_epoch = round(n_training_samples/batch_size)

	n_val_samples = len(x_val)
	val_step_per_epoch = round(n_val_samples/batch_size)

	for epoch in range(epochs):
		print("\nStart of epoch %d" % (epoch,))

		train_loss = []   
		val_loss = []

		training_steps = 0
		# Iterate over the batches of the dataset
		for (x_batch_train, y_batch_train) in generator_data_average_profile(x_train,
																			y_train,
																			chr_seq,
																			input_augmentation,
																			batch_size):
			loss_value = train_step(x_batch_train, y_batch_train, 
									model, optimizer, 
									loss_fn, train_acc_metric)
			train_loss.append(float(loss_value))
			training_steps+=1
			if (training_steps%training_step_per_epoch == 0): break


		# Run a validation loop at the end of each epoch
		val_steps = 0
		for (x_batch_val, y_batch_val) in generator_data_average_profile(x_val,
																			y_val,
																			chr_seq,
																			input_augmentation,
																			batch_size):
			val_loss_value = test_step(x_batch_val, y_batch_val, 
									   model, loss_fn, 
									   val_acc_metric)
			val_loss.append(float(val_loss_value))

			val_steps+=1
			if (val_steps%val_step_per_epoch == 0): break
			
		# Display metrics at the end of each epoch
		train_acc = train_acc_metric.result()
		print("training RMSE over epoch: %.4f" % (float(train_acc),))

		val_acc = val_acc_metric.result()
		print("Validation RMSE : %.4f" % (float(val_acc),))

		# Reset metrics at the end of each epoch
		train_acc_metric.reset_states()
		val_acc_metric.reset_states()

		# log metrics using wandb.log
		wandb.log({'epochs': epoch,
				   'loss': np.mean(train_loss),
				   'RMSE': float(train_acc), 
				   'val_loss': np.mean(val_loss),
				   'val_RMSE':float(val_acc)})