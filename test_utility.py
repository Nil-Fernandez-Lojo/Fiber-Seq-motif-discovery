import numpy as np
from util import generator_data,break_reads_to_fixed_size

seq_pos = np.array([
[0,0,30],
[0,30,60],
[0,90,120],
[0,60,90],
[0,150,180],
[0,120,150],
[0,180,210],
[0,210,240],
[0,0,30]])
n_methyl = np.cumsum(np.array([0,1,2,3,0,1,2,3,0,0]))
methyl_rel_pos = np.array([0,10,20,0,10,20,0,10,20,0,14,19])
chr_seq = ['ACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatanACtTgNnGATccnGcattaccggnTAgnatctatnngattcccatan']
augmentation_value = 2
batch_size = 3
output_length = 18
print('Seq pos')
print(seq_pos)
print('n_methyl')
print(n_methyl)
print('methyl_rel_pos')
print(methyl_rel_pos)
seq_pos,n_methyl,methyl_rel_pos = break_reads_to_fixed_size(seq_pos,n_methyl,methyl_rel_pos,output_length)
print('After fixing size:')
print('Seq pos')
print(seq_pos)
print('n_methyl')
print(n_methyl)
print('methyl_rel_pos')
print(methyl_rel_pos)

# i = 0
# for batch_input, batch_output in generator_data(seq_pos,n_methyl,methyl_rel_pos, chr_seq,augmentation_value,batch_size):
# 	print("Generator iteration ",i)
# 	print(batch_input)
# 	print(batch_output)
# 	print()
# 	i+=1
# 	if i>6: break
