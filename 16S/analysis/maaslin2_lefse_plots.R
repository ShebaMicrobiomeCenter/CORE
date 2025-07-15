source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(Maaslin2)
library(ggplot2)

taxonomy_flag = T
taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-29_merge/res/16S_DB1-29_merged/taxonomy/taxonomy_gg2_source_risk.tsv'

name = 'CORE'
path = '../res_stool_TI_v3_amnonClean/'
taxa_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, name)
na_str = c('no_data','_','NA','unknown', 'other','na','No_followup')
taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str, quote = '', comment.char = '')

# taxa$gastro_dass21_total_stress = as.numeric(taxa$gastro_dass21_total_stress)

multi_file = '../../multiomics/data/CORE_screening_multiomics_map_v8.txt'
multi_map = read.table(multi_file, header = T, sep = '\t', na.strings = c('NA','','#N/A'), comment.char = '',quote = '')

multi_map = multi_map[!is.na(multi_map$ID_16s_stool),]
row.names(multi_map) = multi_map$ID_16s_stool

taxa = taxa[taxa$SampleID %in% multi_map$ID_16s_stool,]
multi_map = multi_map[taxa$SampleID,]
## checked age gender Dx activity vs multiomics_map_v8. no problems
# multi_map$Gender == taxa$Gender

head_path = sprintf('%s/R_res/maaslin2/',path)
dir.create(head_path)



# name = 'Group_age_gender_noCtrl'


# # name = 'CORE_age_gender_bmi_location_activity_socioeconomi'
# name = 'CORE_age_gender_bmi_location_activity_peripherality'
# fixed_cols = c('Age','Gender','BMI','Source','Dx_activity',
#                 'city_peripherality_index_2020_cluster')
#                # 'city_socioeconomi_cluster_2019')
# levels_ref = c("Dx_activity,Control")
# # random_cols = c('pn_ID')
# random_cols = c()
# metadata_cols = c(fixed_cols, random_cols)
# x_var = 'SampleID'; group_var = 'Dx_activity'; taxa_sig_var = 'Dx_activity'

# name = 'CORE_age_gender_bmi_location_activity'
# fixed_cols = c('Age','Gender','BMI','Source','Dx_activity')
# levels_ref = c("Dx_activity,Control")
# # random_cols = c('pn_ID')
# random_cols = c()
# metadata_cols = c(fixed_cols, random_cols)
# x_var = 'SampleID'; group_var = 'Dx_activity'; taxa_sig_var = 'Dx_activity'

# name = 'CORE_age_gender_bmi_activity_TI'
# taxa = taxa[taxa$Source == 'TI',]
# # name = 'CORE_age_gender_bmi_activity_stool'
# # taxa = taxa[taxa$Source == 'stool',]
# fixed_cols = c('Age','Gender','BMI','Dx_activity')
# levels_ref = c("Dx_activity,Control")
# # random_cols = c('pn_ID')
# random_cols = c()
# metadata_cols = c(fixed_cols, random_cols)
# x_var = 'SampleID'; group_var = 'Dx_activity'; taxa_sig_var = 'Dx_activity'

# name = 'CORE_only_stool_age_gender_dass21'
# taxa = taxa[taxa$Cohort == 'CORE',]
# # name = 'CORE_age_gender_bmi_activity_stool'
# taxa = taxa[taxa$Source == 'stool',]
# fixed_cols = c('Age','Gender',
#                'gastro_dass21_total_stress')
# # levels_ref = c("Dx_activity,Control")
# # random_cols = c('pn_ID')
# random_cols = c()
# metadata_cols = c(fixed_cols, random_cols)
# x_var = 'SampleID'; group_var = 'Dx_activity'; taxa_sig_var = 'Dx_activity'

name = 'CORE_stool_age_gender_Dx_activity'
fixed_cols = c('Age','Gender','Dx_activity')
levels_ref = c("Dx_activity,Control")
# random_cols = c('pn_ID')
random_cols = c()
metadata_cols = c(fixed_cols, random_cols)
x_var = 'SampleID'; group_var = 'Dx_activity'; taxa_sig_var = 'Dx_activity'

maas_path = sprintf('%s/%s',head_path, name)

ASV = write_maaslin2(taxa, metadata_cols, maas_path, taxonomy_flag = taxonomy_flag, taxonomy_file = taxonomy_file)
# ASV = write_maaslin2(taxa, metadata_cols, maas_path, taxonomy_flag = taxonomy_flag)

out_path = sprintf('%s/res/',maas_path)

fit_data <- Maaslin2( input_data = sprintf('%s/taxonomy.tsv',maas_path),
                      input_metadata = sprintf('%s/metadata.tsv',maas_path),
                      output = out_path,
                      reference = levels_ref,
                      # transform = "AST",
                      fixed_effects = fixed_cols,
                      random_effects = random_cols)

# if changes to taxonomy, add the asv to the result files
if ( taxonomy_flag )
{
  maaslin_add_asv_to_res(ASV, out_path)
}
# maaslin2_plot_colored_scatters(taxa, x_var = taxa_sig_var, col_var = group_var, head_path, name, out_dir_name = sprintf('my_%s_scatters', taxa_sig_var), fdr = 0.25)


p_heatmap25_res = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', taxa_sig_var = taxa_sig_var, fdr = 0.25)
p_heatmap25 = p_heatmap25_res[[1]]
# p_heatmap1 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', taxa_sig_var = taxa_sig_var, fdr = 0.1)

taxa_m = p_heatmap25_res[[2]]
taxa_m = merge(taxa_m, taxa[,c('SampleID','sample_ID',# 'Source',
                               fixed_cols[fixed_cols!=group_var])], 
               by='SampleID',all.x = T)
taxa_m$Dx_activity = factor(taxa_m$Dx_activity, levels = c('Control','CD_remission','CD_Active'))
taxa_m$id = sprintf('%s',taxa_m$sample_ID)
temp=unique(taxa_m$id)
taxa_m$id = factor(taxa_m$id, levels = c( temp[ order( as.numeric(gsub('_.*','',temp)) ) ] ))

library(viridis)
p_heatmap = ggplot(taxa_m, aes_string(x='id', y='taxa', fill = 'log10_RA', colour = 'log10_RA')) +
  geom_tile() +
  theme_classic() +
  ylab('ASV') + # xlab('Days from accident') +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        # axis.text.x = element_text(angle=90, hjust = 1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_viridis() +
  scale_colour_viridis() +
  # facet_grid(as.formula(sprintf('.~%s+%s','Dx_activity','Source')),
  #            scales = 'free', space = 'free', )
  facet_grid(as.formula(sprintf('.~%s','Dx_activity')),
           scales = 'free', space = 'free', )
# out_fig = sprintf('%s/my_heatmaps/%s_sigTaxa_%s_%s_heatmap%s_fdr%s.tiff', out_path, taxa_sig_var, x_var, group_var, '',0.25)
# ggsave(out_fig,plot = p_heatmap, device = 'tiff', width = 5.2,height =2.5, compression  = 'lzw')
# 
p_heatmap2 = p_heatmap + theme(axis.text.x = element_text(angle=90, hjust = 1, size=5))
# out_fig = sprintf('%s/my_heatmaps/%s_sigTaxa_%s_%s_heatmap%s_fdr%s_withDays.tiff', out_path, taxa_sig_var, x_var, group_var, '',0.25)
# ggsave(out_fig,plot = p_heatmap2, device = 'tiff', width = 8,height =3, compression  = 'lzw')
p_heatmap3 = p_heatmap2 + theme(axis.text.y = element_text(angle=0, hjust = 1, size=5))
out_fig = sprintf('%s/my_heatmaps/%s_sigTaxa_heatmap_fdr%s_full.tiff', out_path, taxa_sig_var,0.25)
ggsave(out_fig,plot = p_heatmap3, device = 'tiff', width = 12,height =12, compression  = 'lzw')


cols = c('#0000ff','#ffa500','#ff0000')

temp = unique(taxa_m[,c('taxa','coef')])
temp= temp[order(temp$coef),]
top_num = 10
# top_tx = temp$taxa[1:top_num]
# top_tx2 = temp$taxa[dim(temp)[1]:( dim(temp)[1]-(top_num-1) )]
all_tx1 <- aggregate(coef ~ taxa, data = temp[temp$coef>0,], max)
all_tx2 <- aggregate(coef ~ taxa, data = temp[temp$coef<0,], min)
all_tx1 = all_tx1[order(all_tx1$coef, decreasing = T),]
top_tx = all_tx1$taxa[1:top_num]
all_tx2 = all_tx2[order(all_tx2$coef),]
top_tx2 = all_tx2$taxa[1:top_num]
tx = c(top_tx,top_tx2)
# tx = 'c__Clostridia.o__Clostridiales_ASV15397'
# # tx = unique(taxa_m$taxa[ abs(taxa_m$coef)>3 ])
# tx = c('g__Faecalibacterium.s__prausnitzii_ASV00199',
#        'c__Clostridia.o__Clostridiales_ASV15397',
#        'g__.Ruminococcus..s__torques_ASV15337')
taxa_m$RA = 10^( taxa_m$log10_RA )
taxa_m$lb = gsub('_.*','',taxa_m$id)
taxa_m$lb[taxa_m$lb == 'NA'] = ''
t=unique(taxa_m[,c('taxa','coef')])
for (i in dim(t)[1]:1)
{
  if (sum(t$taxa[i] == t$taxa) & t$coef[i] != max(t$coef[t$taxa == t$taxa[i]]))
    t=t[-i,]
}
taxa_m$taxa = factor(taxa_m$taxa, levels = t$taxa[order(t$coef)])
p_lefse = ggplot(taxa_m[taxa_m$taxa %in% tx,], aes(x=id, y = RA*10000)) +
  # geom_bar(stat="identity", fill='cyan4') +
  # geom_bar(stat="identity",aes(fill=log(as.numeric(lb)))) +
  geom_bar(stat="identity",aes(fill= Dx_activity )) +
  # geom_bar(stat="identity",aes(fill= GROUP), alpha=0.7 ,color = 'black', size=0.1) + scale_fill_manual(values = cols) + 
  theme_bw() +
  scale_fill_manual(values = rev(c("#ff8080", "#ffc080", "#8d8dff"))) + 
  # ylab('RA') + 
  # xlab('Days from accident') +
  scale_x_discrete(labels = taxa_m$lb, breaks = taxa_m$id ) + 
  scale_y_log10(expand = expansion(mult = c(0, 0.1))) +
  # scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  facet_grid(taxa~Dx_activity, scales = 'free', space = 'free_x')+
  theme(axis.text.x = element_text(angle=90, hjust = 1, size=5),
        axis.text.y = element_text(size=5),
        # strip.background = element_rect(fill='white'),
        panel.grid = element_blank(), strip.text.y.right = element_text(angle = 0)) 
  # facet_grid(as.formula(sprinthttp://172.35.15.242:8799/graphics/plot_zoom_png?width=973&height=279f('%s~%s','taxa','GROUP')),
  #            scales = 'free', space = 'free', )
p_lefse
out_fig = sprintf('%s/my_heatmaps/%s_sigTaxa_lefse_top%s.tiff', out_path, taxa_sig_var,top_num)
# ggsave(out_fig,plot = p_lefse, device = 'tiff', width = 12,height =12, compression  = 'lzw')
ggsave(out_fig,plot = p_lefse, device = 'tiff', width = 12,height =6, compression  = 'lzw')

taxa_m2 = unique( taxa_m[,c('taxa','coef')] )
p = p_bar = ggplot(taxa_m2, aes(x=coef, y = taxa, fill=coef<0)) +
  geom_bar(stat="identity") + theme_classic() + theme(axis.text.y = element_text(size=6))
out_fig = sprintf('%s/my_heatmaps/%s_sigTaxa_bar_fdr%s.tiff', out_path, taxa_sig_var,0.25)
ggsave(out_fig,plot = p, device = 'tiff', width = 7,height =12, compression  = 'lzw')

