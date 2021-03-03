from weblogo import *
seq_folder = 'seq_filter_1_weblogo'
logo_folder = 'weblogo_filters_layer_1'

for i in range(64):
	fin = open(seq_folder+'/filter_'+str(i)+'.fa')
	seqs = read_seq_data(fin)
	logodata = LogoData.from_seqs(seqs)
	logooptions = LogoOptions()
	logooptions.title = "Weblogo filter "+ str(i)
	logoformat = LogoFormat(logodata, logooptions)
	eps = eps_formatter(logodata, logoformat)
	f=open(logo_folder+'/filter_'+str(i)+'.eps', 'wb')
	f.write(eps)
	f.close()