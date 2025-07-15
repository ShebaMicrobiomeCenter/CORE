#!/usr/bin/env python
import os
import pandas as pd 
import sys
sys.path.insert(0, '/pita/users/tzipi/code/rnaSeq')
from run_kallisto_funcs import *

kallisto_index = '/pita/users/tzipi/bin/kallisto_linux-v0.42.5/index/hg_gencode_v24'

# setting variables (paths)
in_path = 'fastq_files/'

#out_path = '/pita/users/tzipi/projects/rnaSeq/celiac_duodenal_plosOne/output/'
out_path = 'kallisto/'
os.system('mkdir ' + out_path)

run_kallisto_paired_kal_42_5_gencode_v24(in_path, out_path, r1_suffix = '_read1.fastq', r2_suffix = '_read2.fastq', kallisto_index = kallisto_index)

'''
# merging transcript level to gene level
gene_out_path = out_path + 'kallisto_gene/'
os.system('mkdir ' + gene_out_path)
trans2gene(kal_out_path, gene_out_path)
'''

