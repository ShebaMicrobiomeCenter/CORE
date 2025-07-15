source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

## clean the original biom file to only "clean" ASVs from Amnon's table,
## while keeping the original read count (non-normalized)
## need this to calculated rarefied alpha diversity for example

## load amnon clean data
name = 'CORE'
path = '../res_stool_TI_v3_amnonClean/'
biom_txt_file = sprintf('%s/biom/feature-table.txt',path)

otu_clean = read.table(file = biom_txt_file, skip = 1, comment.char = '', header = T, row.names = 'OTU_ID')
otu_clean = as.data.frame(t(otu_clean))

## load non-clean and normalized original data
name = 'CORE'
path = '../res_stool_TI_v2/'
biom_txt_file = sprintf('%s/biom/feature-table.txt',path)

otu = read.table(file = biom_txt_file, skip = 1, comment.char = '', header = T, row.names = 'OTU_ID')
otu = as.data.frame(t(otu))

## filter to stool samples
taxa_file = sprintf('%s/biom/%s_deblur_map_L7.1.txt',path, name)
na_str = c('no_data','_','NA','unknown', 'other','na')
taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str, comment.char = '', quote = '')
taxa = taxa[taxa$Source == 'stool',]

# otu_clean = otu_clean[taxa$SampleID,]

## take fitting samples and ASVs
otu = otu[taxa$SampleID, names(otu_clean)]
## remove missing ASVs (existed only in the not-used TI samples)
otu = otu[, as.numeric(colSums(otu))!=0 ]

write.table(t(otu), file = '../data/amnon_sep24/biom_core_stool_by_amnon_clean_ASVs.txt',
            quote = F, sep = '\t', row.names = T, fileEncoding = '')


