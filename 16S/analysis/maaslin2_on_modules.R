source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(Maaslin2)
library(ggplot2)

taxonomy_flag = T
taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-30_merge/res/16S_DB1-30_merged/taxonomy/taxonomy_gg2_source_risk.tsv'
# taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-24_merge/res/16S_DB1-29_merged/taxonomy/metadata_v2.tsv'

name = 'CORE'
path = '../res_stool_TI_v3_amnonClean/'
taxa_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, name)
na_str = c('no_data','_','NA','unknown', 'other','na','No_followup')
taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str, comment.char = '', quote = '')


head_path = sprintf('%s/R_res/maaslin2_modules/',path)
dir.create(head_path)

multi_file = '../../multiomics/data/CORE_screening_multiomics_map_v6.txt'
multi_map = read.table(multi_file, header = T, sep = '\t', na.strings = c('NA','','#N/A'), comment.char = '',quote = '')

taxa = taxa[taxa$SampleID %in% multi_map$ID_16s_stool,]

## merge WGCNA like module data to be used in the maaslin analysis
MEs_file = '/pita/users/tzipi/projects/multiomics/CORE/rnaSeq/res/WGCNA_like_enrichemnt/CORE_S_plus/CORE_S_plus_by_enrichment_genes_module_MEs.txt'
me = read.table(MEs_file, header = T, sep = '\t', na.strings = c('NA','','#N/A'), comment.char = '',quote = '')
temp = multi_map[,c('ID_16s_stool','ID_rnaSeq_TI')]
temp$ID_rnaSeq_TI = make.names(temp$ID_rnaSeq_TI)
temp = temp[temp$ID_rnaSeq_TI %in% row.names(me) &
              !is.na(temp$ID_rnaSeq_TI) & 
              !is.na(temp$ID_16s_stool), ]
taxa = taxa[taxa$SampleID %in% temp$ID_16s_stool, ]
# organize order
row.names(temp) = temp$ID_16s_stool
temp = temp[taxa$SampleID, ]
me = me[temp$ID_rnaSeq_TI, ]

taxa_old = taxa
module_names = names(me)
for (module in module_names)
{
  taxa = taxa_old
  names(taxa)[names(taxa) == 'ID_MBX_serum'] = module
  taxa[[module]] = me[[module]]
  name = sprintf('CORE_stool_age_gender_%s', module)
  fixed_cols = c('Age','Gender',module)
  # name = sprintf('CORE_stool_CD_age_gender_%s', module)
  # taxa = taxa[taxa$Dx_activity != 'Control', ]
  # fixed_cols = c('Age','Gender',module)
  # name = sprintf('CORE_stool_noActive_age_gender_%s', module)
  # taxa = taxa[taxa$Dx_activity != 'CD_Active', ]
  # fixed_cols = c('Age','Gender',module)
  # name = sprintf('CORE_stool_rem_age_gender_%s', module)
  # taxa = taxa[taxa$Dx_activity == 'CD_remission', ]
  # fixed_cols = c('Age','Gender',module)
  # name = sprintf('CORE_stool_age_gender_Dx_activity_%s', module)
  # fixed_cols = c('Age','Gender','Dx_activity',module)
  # levels_ref = c("Dx_activity,Control")
  # random_cols = c('pn_ID')
  random_cols = c()
  metadata_cols = c(fixed_cols, random_cols)
  x_var = 'SampleID'; group_var = 'Gender'; taxa_sig_var = 'Gender'
  
  
  maas_path = sprintf('%s/%s',head_path, name)
  
  ASV = write_maaslin2(taxa, metadata_cols, maas_path, taxonomy_flag = taxonomy_flag, taxonomy_file = taxonomy_file)
  # ASV = write_maaslin2(taxa, metadata_cols, maas_path, taxonomy_flag = taxonomy_flag)
  
  out_path = sprintf('%s/res/',maas_path)
  
  fit_data <- Maaslin2( input_data = sprintf('%s/taxonomy.tsv',maas_path),
                        input_metadata = sprintf('%s/metadata.tsv',maas_path),
                        output = out_path,
                        # reference = levels_ref,
                        # transform = "AST",
                        fixed_effects = fixed_cols,
                        random_effects = random_cols)
  
  # if changes to taxonomy, add the asv to the result files
  if ( taxonomy_flag )
  {
    maaslin_add_asv_to_res(ASV, out_path)
  }
  
  
}

# maaslin2_plot_colored_scatters(taxa, x_var = taxa_sig_var, col_var = group_var, head_path, name, out_dir_name = sprintf('my_%s_scatters', taxa_sig_var), fdr = 0.25)


p_heatmap25_res = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', taxa_sig_var = taxa_sig_var, fdr = 0.25)
p_heatmap25 = p_heatmap25_res[[1]]
# p_heatmap1 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', taxa_sig_var = taxa_sig_var, fdr = 0.1)

# taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str, comment.char = '', quote = '')
# taxa$Group[taxa$Cohort == 'IIRN'] = 'IIRN'
# taxa$group2_cohort = sprintf('%s_%s', taxa$Cohort, taxa$group2)
# taxa = taxa[taxa$group2 %in% c('Never_flared','pre_Flared') & !is.na(taxa$group2),]
# group_var = 'group2_cohort'
# p_heatmap25 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_all_samples', taxa_sig_var = taxa_sig_var, fdr = 0.25)

# taxa_m = p_heatmap25_res[[2]]
# taxa_m = merge(taxa_m, taxa[,c('SampleID','GROUP','BMI','days_from_accident')], by='SampleID',all.x = T)
# taxa_m$GROUP = factor(taxa_m$GROUP, levels = c('Control','SUB_ACUTE','CHRONIC'))
# # taxa_m$id = sprintf('%s_%s',taxa_m$BMI, taxa_m$SampleID)
# taxa_m$id = sprintf('%s_%s',taxa_m$days_from_accident, taxa_m$SampleID)
# temp=unique(taxa_m$id)
# taxa_m$id = factor(taxa_m$id, levels = c( temp[ order( as.numeric(gsub('_.*','',temp)) ) ] ))
# 
# library(viridis)
# p_heatmap = ggplot(taxa_m, aes_string(x='id', y='taxa', fill = 'log10_RA', colour = 'log10_RA')) + 
#   geom_tile() + 
#   theme_classic() + 
#   ylab('ASV') + xlab('Days from accident') + 
#   scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#         # axis.text.x = element_text(angle=90, hjust = 1), 
#         axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
#   scale_fill_viridis() + 
#   scale_colour_viridis() + 
#   facet_grid(as.formula(sprintf('.~%s','GROUP')), 
#              scales = 'free', space = 'free', ) 
# out_fig = sprintf('%s/my_heatmaps/%s_sigTaxa_%s_%s_heatmap%s_fdr%s.tiff', out_path, taxa_sig_var, x_var, group_var, '',0.25)
# ggsave(out_fig,plot = p_heatmap, device = 'tiff', width = 5.2,height =2.5, compression  = 'lzw')
# 
# p_heatmap2 = p_heatmap + theme(axis.text.x = element_text(angle=90, hjust = 1, size=5))
# out_fig = sprintf('%s/my_heatmaps/%s_sigTaxa_%s_%s_heatmap%s_fdr%s_withDays.tiff', out_path, taxa_sig_var, x_var, group_var, '',0.25)
# ggsave(out_fig,plot = p_heatmap2, device = 'tiff', width = 8,height =3, compression  = 'lzw')




# taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)
# p_heatmap25 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_all_samples', taxa_sig_var = taxa_sig_var, fdr = 0.25)
# p_heatmap1 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_all_samples', taxa_sig_var = taxa_sig_var, fdr = 0.1)
# 

# path = '../res_israel_noRS_amnonFlt/'
# taxa_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, 'source')
# taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)
# p_heatmap25 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_israel_noRS', taxa_sig_var = taxa_sig_var, fdr = 0.25)
# p_heatmap1 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_israel_noRS', taxa_sig_var = taxa_sig_var, fdr = 0.1)

# fdr = 0.25
# taxa_m = maaslin2_set_data_for_heatmap_plot(taxa, head_path, name, x_var, taxa_sig_var, group_var = group_var, taxonomy_file = taxonomy_file, fdr)
# taxa_m = merge(taxa_m, taxa[,c('SampleID','location')], by = 'SampleID', all.x=T, sort = F)
# 
# library(viridis)
# p_heatmap = ggplot(taxa_m, aes_string(x=x_var, y='taxa', fill = 'log10_RA')) + geom_tile() + 
#   facet_grid(~Patient_group2+location, scales = 'free') + 
#   theme_classic() + 
#   ylab('ASV') + 
#   scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   scale_fill_viridis() 
# out_path = sprintf('%s/%s/res/%s',head_path,name, 'my_heatmaps')
# dir.create(out_path)
# out_fig = sprintf('%s/%s_sigTaxa_%s_%s_heatmap%s_fdr%s.tiff', out_path, taxa_sig_var, x_var, group_var, '_israel_noRS', fdr)
# ggsave(out_fig,plot = p_heatmap, device = 'tiff', width = 12,height = 2+length(unique(taxa_m$ASV))/8, compression  = 'lzw')

