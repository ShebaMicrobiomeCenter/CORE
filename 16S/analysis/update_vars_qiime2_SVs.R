
library(ggplot2)
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

name = 'CORE'
path = '../res_stool_TI_v3_amnonClean/'
# path = '../res_stool_v3_amnonClean/'

# start with biom table and updates taxonomic taxa file. creates an updated SVs file (with metadata, aDiv etc.)

biom_file = sprintf('%s/biom/feature-table.biom',path)
biom_txt_file = sprintf('%s/biom/feature-table.txt',path)
sv_file = sprintf('%s/biom/%s_deblur_map_SVs.txt',path, name)
out_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, name)
# out_file = sprintf('%s/biom/%s_deblur_map_SVs.1_nonNorm.txt',path, name)

taxa_file = sprintf('%s/biom/%s_deblur_map_L7.1.txt',path, name)
taxa_out_file = sprintf('%s/biom/%s_deblur_map_L7.2.txt',path, name) # with health index

# convert biom to text
# system(sprintf('biom convert -i %s -o %s --to-tsv', biom_file, biom_txt_file))
# convert the biom to dataframe I can work with
system(sprintf('sed -i \'s/#OTU ID/OTU_ID/g\' %s', biom_txt_file))

otu = read.table(file = biom_txt_file, skip = 1, comment.char = '', header = T, row.names = 'OTU_ID')
otu = t(otu)
otu = as.data.frame(otu)
otu$SampleID <- rownames(otu)

#read taxa file
na_str = c('no_data','_','NA','unknown', 'other','na')
taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str, comment.char = '', quote = '')
# names(temp)[1] = 'SampleID'

taxa_old = taxa

# replace taxonomic values with SVs values
# pos = which( names(taxa) == 'Description' )
# metadata_pos = 1:pos
metadata_pos = which( !grepl('k__[AB]',names(taxa)) )
pos = max(metadata_pos)

taxa = taxa[,metadata_pos]
taxa = merge(taxa, otu, by='SampleID')


taxa_pos = (pos+1):dim(taxa)[2]
for ( i in 1:dim(taxa)[1] )
{
  taxa[i,taxa_pos] = taxa[i,taxa_pos]/sum(taxa[i,taxa_pos])
}

## add health index, binary version
hi = calculate_binary_health_index(taxa)
temp = data.frame(SampleID = taxa$SampleID, health_index = hi)
taxa = merge(temp, taxa, by='SampleID',all.x = T)  

## add IBD indeex
ibdi = calculate_health_index(taxa, 
                                     up_file = '/pita/users/tzipi/projects/multiomics/CORE/16S/data/IBD_specific_taxa/down_IBD_ASVs.txt', 
                                     down_file = '/pita/users/tzipi/projects/multiomics/CORE/16S/data/IBD_specific_taxa/up_IBD_ASVs.txt')
temp2 = data.frame(SampleID = taxa$SampleID, IBD_index = ibdi)
taxa = merge(temp2, taxa, by='SampleID',all.x = T)  


write.table( taxa, file = out_file, sep = '\t',quote = F,row.names = F)

## rewriting
taxa_old = merge(temp, taxa_old, by='SampleID',all.x = T)  
taxa_old = merge(temp2, taxa_old, by='SampleID',all.x = T)  
write.table( taxa_old, file = taxa_out_file, sep = '\t',quote = F,row.names = F)

