#!/usr/bin/env python
import os
import glob
import pandas as pd 
import gzip

# fastq_path = '/pita/users/tzipi/data/rnaSeq/RISK_lexogen/SK' + run_num + '_merged/'
# kal_path = './res/kallisto_42.5_gencode_v24/SK' + run_num + '_merged/'

fastq_path = 'fastq_human/'
kal_path = 'kallisto_human/'

## getting number of reads in the fastq files (before mapping)
files = glob.glob(fastq_path + '*read1.fastq*')
names = [ file.split('/')[-1].split('_read1.')[0] for file in files ]
read_num = [0] * len(names)
for i in range(0,len(files)):
	if files[i].split('.')[-1] == 'gz':
		read_num[i] = sum(1 for line in gzip.open(files[i])) /4
	elif files[i].split('.')[-1] == 'fastq':
		read_num[i] = sum(1 for line in open(files[i])) /4
	else:
		print('Invalid file type for ' + files[i])
		read_num[i] = 0
	

kal_read_num = [0] * len(names)
## getting number of reads in the kallisto output
kal_files = [ kal_path + name + '/abundance.tsv' for name in names] 
for i in range(0,len(kal_files)):
	df = pd.read_csv(kal_files[i], sep = '\t')
	kal_read_num[i] = int(sum(df['est_counts']))

res_df = list(zip( names, read_num, kal_read_num ))
res_df = pd.DataFrame(data = res_df, columns=['Sample', 'Reads_fastq','Reads_kallisto_res'])
res_df.to_csv('reads_count.tsv',index=False,header=True, sep = '\t')

