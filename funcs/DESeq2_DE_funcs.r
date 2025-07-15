
library('DESeq2')

# read the htseq files by metadata samples row, and set DESeqDataSet
read_htseq_files = function(path, metadata_file, condition_col = 'diagnosis',id_col = 'Samples')
{
	# sampleFiles = grep('abundance',list.files(path),value=TRUE)
	metadata = read.table(metadata_file,sep="\t", header=TRUE)

	samples = c()
	for ( i in 1:length(metadata[[id_col]] ) )
	{
		# samples[i] = sprintf('%s_abundance.txt',metadata$Samples[i])
		temp = grep(metadata[[id_col]][i], list.files(path), value = TRUE)
		if ( length( temp )!=1 )
			temp=temp[1]
		samples[i] = temp
	}
	sampleFiles = samples

	sampleCondition <- metadata[[condition_col]]
	sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
	 
	ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=path, design=~condition)
	
	return(ddsHTSeq)
}

# write the input data into 1 matrix, normalized or read counts
write_input_data = function(dds, out_path, norm_flag = F)
{
	foo <- counts(dds, normalized = norm_flag)
	
	if (norm_flag)
	{
		out_path = sprintf('%s_norm.csv', out_path)
	} else
	{
		out_path = sprintf('%s.csv', out_path)
	}
	
	write.csv(foo, file=out_path)
	
	return(out_path)
}

filter_FC_result = function( result ,FDR_cutoff = 0.05, FC_cutoff = log2(1.5) )
{
  pos = vector(mode = 'logical',length = nrow(result) )
  for ( i in 1:length(pos) )
  {
    pos[i] = abs(result$log2FoldChange[i]) >= FC_cutoff & result$padj[i] <= FDR_cutoff
  }
  result = result[pos & !is.na(pos),]
  return(result)
}

# write the results to a csv file (with p. values, fold change etc. if plot_MA_falg also plot MA)
write_DE_result = function(dds, name = 'temp', filter_flag = F, FDR_cutoff = 0.05, FC_cutoff = log2(1.5) , out_path = 'csv_files/', type_flag = T)
{
  	res<-results(dds)
  	res<-res[order(res$padj),]
  	
  	mcols(res,use.names=TRUE)
  	t_res = as.data.frame(res)
  	
  	# adding Gene row to df (instead of as a row name) and gene type
  	t_res[['Gene']] = row.names(t_res)
  	row.names(t_res) = NULL
  	if ( type_flag )
  	  t_res = add_gene_type(t_res, gene_col = 'Gene')
  	
  	if ( filter_flag )
  	{
  	  t_res = filter_FC_result(t_res, FDR_cutoff, FC_cutoff)
  	  write.table(t_res,file = sprintf('%s%s_deseq2_res_filtered.tsv',out_path, name),quote = F,sep = '\t',row.names = F)
  	} else 
  	{
  	  write.table(t_res,file = sprintf('%s%s_deseq2_res.tsv',out_path, name), quote = F,sep = '\t',row.names = F) 
  	} 
  	return(t_res)
}

plot_DE_MA = function(dds, name = 'temp')
{
  plotMA(dds,ylim=c(-5,5),main='DESeq2')
  dev.copy(png,sprintf('%s_MAplot.png',name))
  dev.off()
}

# claculate DE from htseq files
calc_DE_DESeq2 = function(directory, metadata_file, name, condition_col, condition_levels)
{
	ddsHTSeq = read_htseq_files(directory, metadata_file, condition_col = condition_col)

	colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=condition_levels)

	dds<-DESeq(ddsHTSeq)
	res<-results(dds)
	res<-res[order(res$padj),]

	data_path =  sprintf('$s_data',name)
	write_input_data(dds,data_path, norm_flag = F)

	write_DE_result(dds, name, T)
	
	return(dds)
}

filter_count_matrix_samples = function(count_matrix, samples_list, startX_flag = T)
{
  temp = as.character(samples_list)
  for (i in 1:length(temp) ) 
  {
    if (startX_flag) 
      temp[i] = sprintf('X%s',temp[i]) # if the name starts with a number this is how R reads it.
    # replace '-' with '.'
    temp[i] = gsub(temp[i],pattern = '-',replacement = '.',fixed = T)
  }
  # browser()
  count_matrix <- subset(count_matrix, select = temp) #bad!!!
  return( count_matrix)
}

filter_count_matrix_genes = function(count_matix, TPM_matix, TPM_val_cutoff, TPM_per_cutoff)
{
  pos = vector(mode = 'logical',length = nrow(count_matix) )
  for ( i in 1:nrow(count_matix) )
  {
    temp = as.numeric(TPM_matix[i,])
    pos[i] = sum(temp > TPM_val_cutoff)/length(temp) > TPM_per_cutoff 
  }
  count_matix = count_matix[pos,]
  return(count_matix)
}

# add_gene_type = function(df,  ref_file = 'data/gencode_v23_tran_details.txt', gene_col = 'Gene', gene_pos=6, type_pos=8)
add_gene_type = function(df,  ref_file = 'data/biomart_gene_type.txt', gene_col = 'Gene', gene_pos=1, type_pos=3)
# add_gene_type = function(df,  ref_file = 'data/biomart_HS_genes_GRCh38.p5.txt', gene_col = 'Gene', gene_pos=1, type_pos=3)
{
  # ref = read.table(ref_file,sep="\t", header=FALSE)
  ref = read.table(ref_file,sep="\t", header=TRUE)
  gene = as.character(ref[,gene_pos])
  g_type = ref[,type_pos]
  
  d_genes = df[[gene_col]]
  
  d_type = vector(mode = 'character',length = length(d_genes))
  for ( i in 1:length(d_genes) )
  {
    g = d_genes[i]
    pos =  which( as.character(g)==gene )[1]
    d_type[i] = as.character(g_type[pos])
  }
  df[['gene_type']] = d_type
  return(df)
}

# ## parameteres
# directory = '/nadata/users/tzipi/mid_files/rnaSeq_RISK/kallisto_gencodeV23_gene_htseqLike/'
# # metadata_file = '/home/tzipi/code/rnaSeq_DE/RISK_rnaSeq_metadata.txt'
# metadata_file = '/home/tzipi/code/rnaSeq_DE/RISK_rnaSeq_metadata_nocCD.txt'
# # metadata = read.table(metadata_file,sep="\t", header=TRUE)
# condition_col = 'diagnosis'
# condition_levels = c('Y_Not_IBD','Y_iCD')
# 
# name = 'RISK_diagnosis'



# dds = calc_DE_DESeq2 = function(directory, metadata_file, name, condition_col, condition_levels)


