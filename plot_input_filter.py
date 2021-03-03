import numpy as np
import tensorflow as tf
import seaborn as sns
import matplotlib.pylab as plt

model = tf.keras.models.load_model('model_trained',compile=False)
filter_values = model.layers[1].get_weights()[0]
vmin = np.amin(filter_values)
vmax = np.amax(filter_values)
#vmin = min(np.amin(filter_values), -np.amax(filter_values))
#vmax = max(np.amax(filter_values), -np.amin(filter_values))

if np.abs(vmin)<np.abs(vmax):
	filter_values[filter_values<0] *= np.abs(vmax/vmin)
	max_val = vmax
else:
	filter_values[filter_values>0] *= np.abs(vmin/vmax)
	max_val = -vmin

#filter_values.shape[2]
fig, axes = plt.subplots(8, 8, figsize=(18, 10),sharex=True,sharey=True)
cbar_ax = fig.add_axes([.91, .3, .03, .4])

cmap = sns.color_palette("coolwarm", as_cmap=True)
for i, ax in enumerate(axes.flat):
    sns.heatmap(data=np.transpose(filter_values[:,:,i]),
    			ax=ax,
                cbar=i == 0,
                yticklabels = ['A','T','G','C'],
                vmin = -max_val,
				vmax = max_val,
				cmap = cmap,
                cbar_ax=None if i else cbar_ax)
cbar_ax.set_yticklabels(["-0.4", "-0.2", "0", str(round(0.2*np.abs(vmax/vmin),2)),str(round(0.4*np.abs(vmax/vmin), 2))])

plt.show()