library(tximport)

# library('DESeq2')
source('/pita/users/tzipi/code/rnaSeq/DESeq2_DE_funcs.r')

Group_prm = 'Dx_activity'
two_groups = c('Control','CD_Active')
# two_groups = c('Control','CD_remission')
# two_groups = c('CD_remission','CD_Active')

# Group_prm = 'CORE_outlier_compare'
# two_groups = c('CORE','CORE_outliers')

name = 'CORE_S_plus' 

metadata_file = 'data/CORE_screening_rnaSeq_map.txt' 
df = read.table(metadata_file, header = T, sep = '\t', comment.char = '',quote = '')

# make a transcript to gene invertion table based on a random kallisto-output file
# example_file = sprintf('%s/%s/abundance.tsv', input_dir, files_names[1])
example_file = sprintf('%s/%s/abundance.tsv', df$path[1], df$SampleID[1])
exp = read.table(example_file,sep="\t", header=T)
trans_names = exp$target_id
gene_name = as.character(trans_names)
ens_name = gene_name
for ( i in 1:length(gene_name) )  {   gene_name[i] = strsplit(gene_name[i], '|', fixed = TRUE)[[1]][6] }
for ( i in 1:length(ens_name) )  {   ens_name[i] = strsplit(ens_name[i], '|', fixed = TRUE)[[1]][2] }
GENEID = sprintf('%s__%s', gene_name, ens_name)
# gene_name = gsub(pattern = '_.*',replacement = '',x = gene_name)
tx2gene = data.frame(TXNAME = trans_names, GENEID = GENEID)



# df = data.frame(SampleID = files_names, files=
#                   sprintf('%s/%s/abundance.h5',input_dir, files_names))

df$files = sprintf('%s/%s',df$path, df$SampleID)
df$files = sprintf('%s/abundance.h5',df$files)
df = df[order(df$files),]

# run tximport to summarize transcript to gene level
txi <- tximport(df$files, type = "kallisto", tx2gene = tx2gene)

temp = as.data.frame(txi$abundance); names(temp) = df$SampleID
temp = data.frame(Gene = row.names(temp),temp)
# write.table(temp, file = sprintf('res/%s_txi_res.txt',name),quote = F, sep='\t',row.names = F)

# temp = as.data.frame(txi$counts); names(temp) = df$new_name
# temp = data.frame(Gene = row.names(temp),temp)
# write.table(temp, file = sprintf('res/%s_txi_counts_res_v2.txt',name),quote = F, sep='\t',row.names = F)

# filtering tximport results to wanted genes (protein coding from the table yael send, TPM1>0.2% samples)
gene_file = '/pita/users/tzipi/projects/rnaSeq/lncRNA/tables/from_yael/19815_Protein_coding_genes_GencodeV24.txt'
gene_df = read.table(gene_file, header = T)
# wanted_genes = gene_df$Gene
wanted_genes = sprintf('%s__%s', gene_df$Gene, gene_df$Ensembl_Gene_ID)
txi_genes = row.names(txi$abundance)
wanted_genes_pos1 = txi_genes %in% wanted_genes


## filter TPM to 1 TPM>0.2%
tpm_genes = row.names(txi$abundance)
wanted_genes_pos2 = vector(mode = 'logical',length = length(tpm_genes))
for ( i in 1:length(tpm_genes) )
{
  per = sum(txi$abundance[i,] > 1) / dim(txi$abundance)[2]
  wanted_genes_pos2[i] = per > 0.2
}
wanted_genes_pos = wanted_genes_pos1 & wanted_genes_pos2

txi$abundance = txi$abundance[wanted_genes_pos, ]
txi$counts = txi$counts[wanted_genes_pos, ]
txi$length = txi$length[wanted_genes_pos, ]
# 
name = sprintf('%s_geneFiltered', name)
temp = as.data.frame(txi$abundance); names(temp) = df$SampleID
temp = data.frame(Gene = row.names(temp),temp)
# write.table(temp, file = sprintf('res/%s_txi_res_geneFiltered.txt',name),quote = F, sep='\t', row.names = F)


name = sprintf('%s_%s_(n=%s)_vs_%s_(n=%s)',name, two_groups[2],
               sum(df[[Group_prm]] == two_groups[2]),two_groups[1], sum(df[[Group_prm]] == two_groups[1]))
pos = df[[Group_prm]] %in% two_groups

df = df[pos,]
txi$abundance = txi$abundance[, pos]
txi$counts = txi$counts[, pos]
txi$length = txi$length[, pos]

coldata = data.frame(row.names=df$new_name, diagnosis=as.factor(df[[Group_prm]]) )


ddsHTSeq <- DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~ diagnosis)
colData(ddsHTSeq)$diagnosis<-factor(colData(ddsHTSeq)$diagnosis, levels=c(two_groups[1],two_groups[2]))


dds<-DESeq(ddsHTSeq)

res<-results(dds)
res<-res[order(res$padj),]
head(res)
#
# dir.create('res')
write_DE_result(dds, name = name, out_path = 'res/', type_flag = F)
write_DE_result(dds, name = name, filter_flag = T, FDR_cutoff = 0.05, FC_cutoff = log2(1.5), out_path = 'res/', type_flag = F)


