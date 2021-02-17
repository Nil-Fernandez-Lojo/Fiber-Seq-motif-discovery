import numpy as np

base_to_int_map = {'a' : 0, 'A' : 0, 't': 1, 'T':1, 'g' : 2, 'G' : 2, 'c' : 3, 'C' : 3, 'n' : -1, 'N': -1}

#### Functions to read the data
def load_genome(fasta_file):
	f = open(fasta_file, "r")
	fl = f.readlines()
	f.close()
	sequences = []
	sequences_name = []
	chr_number = 0
	for line in fl:
		if line[0] == ">":
			sequences_name.append(line[1:-1])
			if (chr_number!=0): sequences.append(my_seq)
			my_seq = ''
			chr_number+=1
		else:
			my_seq+=line[:-1]
	sequences.append(my_seq)
	return sequences_name,sequences

def load_bed_file(bed_file,chr_names):
	#For now cannot deal with headers, chr_names must be a list of of strings where each entry is the name of a chromosome

	#Count number of lines	
	f = open(bed_file, "r")
	n_lines = sum(1 for i in f)
	f.close()

	#get input sequences NN	(genomic position)
	f = open(bed_file, "r")
	input_dat = np.zeros((n_lines,3), dtype=int)
	n_methyl = np.zeros(n_lines+1, dtype=int) #cumulative number of methylations per read
	counter_methyl = 0
	for i in range(n_lines):
		line = f.readline()
		line_split = line.strip().split()
		counter_methyl += int(line_split[9])
		input_dat[i,0] = chr_names.index(line_split[0])
		input_dat[i,1] = int(line_split[1]) #start
		input_dat[i,2] = int(line_split[2]) #end
		n_methyl[i+1] = counter_methyl
	f.close()

	#get output sequence NN	(methylation position)
	f = open(bed_file, "r")
	methyl_rel_pos = np.zeros(counter_methyl, dtype=int)
	for i in range(n_lines):
	 	line = f.readline()
	 	line_split = line.strip().split()
	 	methyl_rel_pos[n_methyl[i]:n_methyl[i+1]] = np.array(list(map(int,line_split[11].split(","))),dtype=int)
	f.close()
	return input_dat,n_methyl,methyl_rel_pos

#### Process data for NN
def seq_to_1_hot(seq):
	# This is not vectorised!!!!	
	idxs = [base_to_int_map[b] for b in seq]
	x = np.zeros((len(seq),4))
	for i in range(len(seq)):
		if (idxs[i]!=-1): x[i,idxs[i]] = 1
	return x

def seq_pos_to_1_hot(seq, start, end):
	#For now, if the base is not known, all zeros. 
	#If the segments goes after the end or before the beginning of the chr after augmentation, 
	#those "virtual bases" are also set to all zeros
	subseq = seq[max(start,0):min(end,len(seq))]
	x = seq_to_1_hot(subseq)
	if (start<0): x = np.concatenate((np.zeros((-start,4),dtype=int), x), axis=0)
	if (end>len(seq)-1): x = np.concatenate((x,np.zeros((end-len(seq),4),dtype=int)), axis=0)
	return x

def generator_data(seq_pos,n_methyl,methyl_rel_pos, chr_seq,augmentation_value,batch_size):
	# seq_pos : N x 3 numpy array: where N is the number of reads.
	#			First column indicates the chromosome number in the chr_seq array
	#			Second column indicates the starting position of the read
	#			Third column indicates the ending position of the read 
	#			The sequence of the read i is therefore given by chr_seq[seq_pos[i,0]][seq_pos[i,1]: seq_pos[i,2]]
	# n_methyl: numpy array of length N+1: Cumulative number of methiylations per read
	#			It is used as a pointer for the methyl_rel_pos array
	#			The list of reltive positions of methylations of read i is given from position n_methyl[i] to n_methyl[i+1] in methyl_rel_pos
	# methyl_rel_pos: 1D numpy array with all the relative 
	# chr_seq: List of chromosome sequence given as a string.
	# augmentation_value: int : number of bases added to the input sequences of the both left and right 
	# batch_size: int: number of samples given at each iteration of this generator function


	output_length = seq_pos[0,2] - seq_pos[0,1]
	input_length = output_length+2*augmentation_value
	batch_input = np.zeros((batch_size, input_length, 4))
	batch_output = np.zeros((batch_size,output_length))
	i = 0
	while True:
		for read_number in range(len(seq_pos)): #The reads should have been shuffled prior to this operation!
			chr_indx = seq_pos[read_number,0]
			read_s = seq_pos[read_number,1]
			read_e = seq_pos[read_number,2]
			batch_input[i] = seq_pos_to_1_hot(chr_seq[chr_indx],read_s-augmentation_value, read_e + augmentation_value)
			batch_output[i][methyl_rel_pos[n_methyl[read_number]:n_methyl[read_number+1]]] = 1

			i+=1

			if(i%batch_size == 0):
				yield (batch_input, batch_output)
				i=0
				batch_input = np.zeros((batch_size, input_length, 4))
				batch_output = np.zeros((batch_size,output_length))

def break_reads_to_fixed_size(seq_pos,n_methyl,methyl_rel_pos,fixed_size):
	n_samples = np.sum(np.floor((seq_pos[:,2]-seq_pos[:,1])/fixed_size).astype('int32'))
	seq_pos_new = np.zeros((n_samples,3), dtype=int)
	n_methyl_new =  np.zeros(n_samples+1, dtype=int) 
	methyl_rel_pos_new = np.zeros(len(methyl_rel_pos), dtype=int)

	counter_input = 0
	tot_methyls = 0
	for read_number in range(len(seq_pos)):
		n_methyls_observed_read = 0
		counter_input_start_read = counter_input
		for i in range((np.floor((seq_pos[read_number,2]-seq_pos[read_number,1])/fixed_size)).astype('int32')):

			seq_pos_new[counter_input,0] = seq_pos[read_number,0]
			seq_pos_new[counter_input,1] = seq_pos[read_number,1] + i*fixed_size
			seq_pos_new[counter_input,2] = seq_pos[read_number,1] + (i+1)*fixed_size
			#print(n_methyl[read_number]+n_methyls_observed, n_methyl[read_number+1])
			for j in range(n_methyl[read_number]+n_methyls_observed_read, n_methyl[read_number+1]):
				if (methyl_rel_pos[j] >= ((i+1)*fixed_size)): break
				methyl_rel_pos_new[tot_methyls] = methyl_rel_pos[j] - i*fixed_size
				n_methyls_observed_read += 1
				tot_methyls+=1

			n_methyl_new[counter_input+1] = tot_methyls

			counter_input+=1
	
	methyl_rel_pos_new = methyl_rel_pos_new[0:tot_methyls]
	return seq_pos_new,n_methyl_new,methyl_rel_pos_new

def shuffle_data(seq_pos,n_methyl,methyl_rel_pos):
	permutation = np.random.permutation(len(seq_pos))

	seq_pos_new = seq_pos[permutation,:]
	n_methyl_new = np.zeros(len(n_methyl), dtype=int)
	methyl_rel_pos_new = np.zeros(len(methyl_rel_pos), dtype=int)
	
	for i in range(len(seq_pos)):
		n_methyl_new[i+1] = n_methyl_new[i]+n_methyl[permutation[i]+1]- n_methyl[permutation[i]]
		methyl_rel_pos_new[n_methyl_new[i]:n_methyl_new[i+1]] = methyl_rel_pos[n_methyl[permutation[i]]:n_methyl[permutation[i]+1]] 

	return seq_pos_new,n_methyl_new,methyl_rel_pos_new

def split_training_testing(seq_pos,n_methyl,methyl_rel_pos,fraction_training):
	n_training = round(fraction_training*len(seq_pos))
	
	seq_pos_training = seq_pos[:n_training,:]
	seq_pos_testing = seq_pos[n_training:,:]
	
	n_methyl_training = n_methyl[:n_training+1]
	n_methyl_testing = n_methyl[n_training:]-n_methyl[n_training]

	methyl_rel_pos_training = methyl_rel_pos[:n_methyl[n_training]]
	methyl_rel_pos_testing = methyl_rel_pos[n_methyl[n_training]:]

	return seq_pos_training,n_methyl_training,methyl_rel_pos_training, seq_pos_testing,n_methyl_testing,methyl_rel_pos_testing

