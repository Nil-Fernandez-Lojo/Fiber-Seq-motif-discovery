import numpy as np
from util import generator_data,break_reads_to_fixed_size,shuffle_data,split_training_testing,generator_data, load_genome,load_bed_file,get_average_profile
import matplotlib.pyplot as plt

genome_file = "dm6.fa"
bed_file = "GSM4411218_tracks_m6A_DS75167.dm6.bed"
bed_file_peaks = 'GSM4411217_DS68346_hotspot_fdr05.bed'
chr_names,chr_seq = load_genome(genome_file)
peaks_pos = load_bed_file_peaks(bed_file_peaks,chr_names,chr_seq)
seq_pos,n_methyl,methyl_rel_pos = load_bed_file(bed_file,chr_names,chr_seq)
methyl_profile, coverage = get_average_profile(seq_pos,n_methyl,methyl_rel_pos, chr_seq)



list_reads_no_methylation = []
list_reads_only_T_methylated = []
list_reads_only_A_methylated = []

for i in range(int(len(seq_pos)/2)):
	if (n_methyl[2*i+1] ==  n_methyl[2*i]) and (n_methyl[2*(i+1)] ==  n_methyl[2*i+1]): list_reads_no_methylation.append(i)
	elif (n_methyl[2*i+1] ==  n_methyl[2*i]): list_reads_only_T_methylated.append(i)
	elif (n_methyl[2*(i+1)] ==  n_methyl[2*i+1]): list_reads_only_A_methylated.append(i)


length_read = []
for read in list_reads_only_A_methylated:
	length_read.append(seq_pos[2*i,2]-seq_pos[2*i,1])
for read in list_reads_only_T_methylated:
	length_read.append(seq_pos[2*i,2]-seq_pos[2*i,1])
plt.hist(length_read)
plt.xlabel("length")


chr_idx = seq_pos[1,0]
methyl_pos_chr = seq_pos[1,1] + methyl_rel_pos[n_methyl[1]:n_methyl[2]]
for base_idx in methyl_pos_chr:
	print(chr_seq[chr_idx][base_idx], end='')
print()

a = []
for i in range(1,len(seq_pos)):
	if (n_methyl[i] == n_methyl[i-1]):
		a.append(seq_pos[i,2] - seq_pos[i,1])
a = np.array(a)
plt.hist(a, bins=20)  # density=False would make counts
plt.xlabel('Read length')
plt.show()

# seq_pos = np.array([
# [0,0,30],
# [0,30,60],
# [0,90,120],
# [0,60,90],
# [0,150,180],
# [0,120,150],
# [0,180,210],
# [0,210,240],
# [0,0,30]])
# n_methyl = np.cumsum(np.array([0,1,2,3,0,1,2,3,0,0]))
# methyl_rel_pos = np.array([0,10,20,0,10,20,0,10,20,0,14,19])
# chr_seq = ['ACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatan']
# augmentation_value = 2
# batch_size = 3
# output_length = 18
# fraction_training = 1/3

# # print('Seq pos')
# # print(seq_pos)
# # print('n_methyl')
# # print(n_methyl)
# # print('methyl_rel_pos')
# # print(methyl_rel_pos)
# seq_pos,n_methyl,methyl_rel_pos = break_reads_to_fixed_size(seq_pos,n_methyl,methyl_rel_pos,output_length)
# print('After fixing size:')
# print('Seq pos')
# print(seq_pos)
# print('n_methyl')
# print(n_methyl)
# print('methyl_rel_pos')
# print(methyl_rel_pos)

# seq_pos_training,n_methyl_training,methyl_rel_pos_training, seq_pos_testing,n_methyl_testing,methyl_rel_pos_testing = split_training_testing(seq_pos,n_methyl,methyl_rel_pos,fraction_training)
# print('training data:')
# print('Seq pos')
# print(seq_pos_training)
# print('n_methyl')
# print(n_methyl_training)
# print('methyl_rel_pos')
# print(methyl_rel_pos_training)

# print('testing data:')
# print('Seq pos')
# print(seq_pos_testing)
# print('n_methyl')
# print(n_methyl_testing)
# print('methyl_rel_pos')
# print(methyl_rel_pos_testing)
# i = 0
# for batch_input, batch_output in generator_data(seq_pos,n_methyl,methyl_rel_pos, chr_seq,augmentation_value,batch_size):
# 	print("Generator iteration ",i)
# 	print(batch_input)
# 	print(batch_output)
# 	print()
# 	i+=1
# 	if i>6: break


# for i in range(len(seq_pos)):
#     for j in range(n_methyl[i], n_methyl[i+1]):
# ...             if (chr_seq[seq_pos[i,0]][seq_pos[i,1] + methyl_rel_pos[j]] not in "AaTt"): 
# ...                     counter_reads_non_AT_methylation +=1
# ...                     break