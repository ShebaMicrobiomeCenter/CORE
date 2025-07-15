source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(Maaslin2)
library(ggplot2)

taxonomy_flag = F
# taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/taxonomy/metadata_v2_source.tsv'

metadata_prms = c('Cohort','Dx','Dx_activity','new_diagnosis',
                  'Age','Gender','BMI','Smoking_YN',
                  'pn_ID')
col_parm = 'Dx_activity'

ffq_file = '../data/CORE_FFQ_with_metadata.txt'
ffq = read.table(ffq_file, header = T, sep = '\t')

taxa = ffq

## filtering to 1 sample per subject
## changing the last col before the ffq data (first_sample in this case)
## to "Description" for technincal reasons
taxa$Smoking_YN = taxa$pn_ID
names(taxa)[names(taxa) == 'Smoking_YN'] = 'Description'

# source('/pita/users/tzipi/code/R_figs/WGCNA_funcs.R')
# taxa = clean_metadata_eampty_noVar(taxa)

head_path = sprintf('../res/maaslin2/')
dir.create(head_path)

names(taxa)[1] = 'SampleID'

name = 'Dx_activity_age_gender_ctlRef'
fixed_cols = c( 'Dx_activity' ,'Age','Gender')
# random_cols = c('pn_ID')
levels_ref = c('Dx_activity,Control')
random_cols = c()
metadata_cols = c(fixed_cols, random_cols)
x_var = 'SampleID'; group_var = 'Dx_activity'; taxa_sig_var = 'Dx_activity'

# name = 'Dx_activity_age_gender_remissionRef'
# fixed_cols = c( 'Dx_activity' ,'Age','Gender')
# # random_cols = c('pn_ID')
# levels_ref = c('Dx_activity,CD_remission')
# random_cols = c()
# metadata_cols = c(fixed_cols, random_cols)
# x_var = 'SampleID'; group_var = 'Dx_activity'; taxa_sig_var = 'Dx_activity'


maas_path = sprintf('%s/%s',head_path, name)

ASV = write_maaslin2(taxa, metadata_cols, maas_path, taxonomy_flag = taxonomy_flag, taxonomy_file = taxonomy_file)
# ASV = write_maaslin2(taxa, metadata_cols, maas_path, taxonomy_flag = taxonomy_flag)

out_path = sprintf('%s/res/',maas_path)

fit_data <- Maaslin2( input_data = sprintf('%s/taxonomy.tsv',maas_path),
                      input_metadata = sprintf('%s/metadata.tsv',maas_path),
                      output = out_path,
                      reference = levels_ref,
                      # transform = "AST",
                      # random_effects = random_cols,
                      fixed_effects = fixed_cols)


# if changes to taxonomy, add the asv to the result files
if ( taxonomy_flag )
{
  maaslin_add_asv_to_res(ASV, out_path)
}

wanted_features = unique(fit_data$results$feature[fit_data$results$metadata == taxa_sig_var & 
                     fit_data$results$qval <= 0.25] )
df = fit_data$results
df = df[df$feature %in% wanted_features, ]
df = df[df$metadata == taxa_sig_var, ]
df$group = sprintf('%s / Control', df$value)
# temp = df[df$value == 'CD_remission', ]
temp = df[df$value == 'CD_Active', ]
temp = temp[order(temp$coef), ]
df$feature = factor(df$feature, levels = temp$feature)

pos = df$coef; pos[df$coef<0] = pos[df$coef<0]-0.18

cols = c( '#e31a4c', '#fbba3a', '#4a65c5')
g_bar = ggplot(df, aes(x=coef, y=feature, fill = group)) + 
  geom_bar(position="dodge", stat="identity", width = 0.8)+
  theme_bw() + theme(panel.grid = element_blank()) +
  scale_fill_manual(values = cols, name = '') + 
  geom_text(
    aes(label = ifelse(qval<=0.25, '*',''), x = pos), # Position text slightly to the right of the bar
    position = position_dodge(width = 0.8),        # MUST match geom_bar's dodge width
    hjust = 0,                                     # Horizontal justification: 0=left, 0.5=center, 1=right
    # We want the text to start after the bar
    vjust = 0.7,                                   # Vertical justification (0.5 = center)
    size =3,                                    # Adjust size as needed
    color = "black"                                # Optional: set color
  ) 
g_bar
out_fig = sprintf('%s/bar_%s.tiff', out_path, taxa_sig_var)
ggsave(out_fig,plot = g_bar, device = 'tiff', width = 6,height = 4.8, compression  = 'lzw')


pos2 = df$coef; pos2[df$coef<0] = pos2[df$coef<0]-0.2

cols = c( '#e31a4c', '#fbba3a', '#4a65c5')
g_bar2 = ggplot(df, aes(y=coef, x=feature, fill = group)) + 
  geom_bar(position="dodge", stat="identity", width = 0.8)+
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = cols, name = '') + 
  geom_text(
    aes(label = ifelse(qval<=0.25, '*',''), y = pos2), # Position text slightly to the right of the bar
    position = position_dodge(width = 0.8),        # MUST match geom_bar's dodge width
    vjust = 0.4,                                     # Horizontal justification: 0=left, 0.5=center, 1=right
    # We want the text to start after the bar
    hjust = 0.5,                                   # Vertical justification (0.5 = center)
    size =3,                                    # Adjust size as needed
    color = "black"                                # Optional: set color
  ) 
g_bar2
out_fig2 = sprintf('%s/bar_%s_flip.tiff', out_path, taxa_sig_var)
ggsave(out_fig2,plot = g_bar2, device = 'tiff', width = 6,height = 4.3, compression  = 'lzw')


# 
# ## plot lefse like 
# wanted_features = unique(fit_data$results$feature[fit_data$results$metadata == taxa_sig_var & 
#                                                     fit_data$results$qval <= 0.25] )
# vrs = names(taxa)[!names(taxa) %in% wanted_features]
# vrs = c( metadata_cols, 'SampleID')
# taxa_m = reshape::melt(data = taxa, id.vars = vrs,measure.vars = wanted_features, variable_name = 'FFQ' )
# 
# temp = unique(fit_data$results[fit_data$results$qval <= 0.25 & 
#                                  fit_data$results$metadata == 'Dx_activity',
#                                c('feature','coef','value')])
# temp= temp[order(temp$coef),]
# 
# temp2 = unique(fit_data$results[fit_data$results$metadata == 'Dx_activity',
#                                 c('feature','coef','value')])
# temp2 = temp2[temp2$feature %in% temp$feature & temp2$value == 'Control',]
# # temp2 = temp2[temp2$feature %in% temp$feature & temp2$value == 'CD_Active',]
# temp2= temp2[rev(order(temp2$coef)),]
# 
# taxa_m$value[taxa_m$FFQ == 'Arachidonic.Cholesterol_mg.d_cal.Choline.Eggs_zmean'] = 
#   taxa_m$value[taxa_m$FFQ == 'Arachidonic.Cholesterol_mg.d_cal.Choline.Eggs_zmean'] + 1.6
# 
# taxa_m$FFQ = factor(taxa_m$FFQ, levels = temp2$feature)
# 
# top_num = 5
# top_tx = temp$feature[1:top_num]
# top_tx2 = temp$feature[dim(temp)[1]:( dim(temp)[1]-(top_num-1) )]
# tx = c(top_tx,top_tx2)
# 
# cols = c( '#e31a4c', '#fbba3a', '#4a65c5')
# # p_lefse = ggplot(taxa_m[taxa_m$FFQ %in% tx,], aes(x=SampleID, y =value)) +
# p_lefse = ggplot(taxa_m, aes(x=SampleID, y =value)) +
#   # geom_bar(stat="identity", fill='cyan4') +
#   # geom_bar(stat="identity",aes(fill=log(as.numeric(lb)))) +
#   # geom_bar(stat="identity",aes(fill= coef<0 )) +
#   scale_fill_manual(values =cols) + 
#   geom_bar(stat="identity", aes(fill =Dx_activity )) +
#   # geom_bar(stat="identity",aes(fill= GROUP), alpha=0.7 ,color = 'black', size=0.1) + scale_fill_manual(values = cols) + 
#   theme_bw() +
#   # ylab('RA') + 
#   xlab('') +
#   # scale_x_discrete(labels = taxa_m$lb, breaks = taxa_m$id ) + 
#   # scale_y_log10(expand = expansion(mult = c(0, 0.1))) +
#   # scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   facet_grid(FFQ~Dx_activity, scales = 'free', space = 'free_x')+
#   theme(# axis.text.x = element_text(angle=90, hjust = 1, size=5),
#         axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#         axis.text.y = element_text(size=5),
#         # strip.background = element_rect(fill='white'),
#         panel.grid = element_blank(), 
#         strip.text.y.right = element_text(angle = 0)) 
# # facet_grid(as.formula(sprinthttp://172.35.15.242:8799/graphics/plot_zoom_png?width=973&height=279f('%s~%s','taxa','GROUP')),
# #            scales = 'free', space = 'free', )
# p_lefse
# # out_fig = sprintf('%s/lefse_like_%s_top10.tiff', out_path, taxa_sig_var)
# # ggsave(out_fig,plot = p_lefse, device = 'tiff', width = 11,height = 3.5, compression  = 'lzw')
# 
# out_fig = sprintf('%s/lefse_like_%s.tiff', out_path, taxa_sig_var)
# ggsave(out_fig,plot = p_lefse, device = 'tiff', width = 11,height = 7, compression  = 'lzw')
# 
# 
# ffqs = unique(taxa_m$FFQ)
# taxa_m$ffq_z = NA
# for ( f in ffqs )
# {
#   pos = which(taxa_m$FFQ == f)
#   taxa_m$ffq_z[pos] = scale(taxa_m$value[pos], center = T, scale = T)
# }
# 
# # ggplot(taxa_m, aes(x = ffq_z, y = FFQ, 
# library(ggridges)
# # p_ridge = ggplot(taxa_m[taxa_m$FFQ %in% tx,], 
# # p_ridge = ggplot(taxa_m[taxa$Gender == 'Female',], 
# p_ridge = ggplot(taxa_m, 
#                  aes(x = ffq_z, y = FFQ, 
#                    color = Dx_activity, 
#                    point_color = Dx_activity, fill = Dx_activity)) +
#   ggridges::geom_density_ridges(
#     # jittered_points = TRUE, 
#     scale = .95, 
#     # rel_min_height = .01, 
#     alpha = 0.7
#     # point_shape = "|", point_size = 3, size = 0.25
#     # position = position_points_jitter(height = 0)
#   ) +
#   scale_y_discrete(expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_fill_manual(values = cols) +
#   scale_color_manual(values = cols, guide = "none") +
#   facet_grid(.~Gender, scale = 'free') + 
#   scale_discrete_manual("point_color", values = cols, guide = "none") +
#   # coord_cartesian(clip = "off") + 
#   theme_ridges(grid = FALSE) 
#   # guides(fill = guide_legend(
#   #   override.aes = list(
#   #     fill = c(cols))
#   # )
#   # ) 
# p_ridge
# # out_fig = sprintf('%s/maaslin_ridges_%s_top10.tiff', out_path, taxa_sig_var)
# # ggsave(out_fig,plot = p_ridge, device = 'tiff', width = 10,height = 5, compression  = 'lzw', bg = "white")
# 
# # out_fig = sprintf('%s/maaslin_ridges_%s.tiff', out_path, taxa_sig_var)
# # ggsave(out_fig,plot = p_ridge, device = 'tiff', width = 10,height = 10, compression  = 'lzw', bg = "white")
# 
# out_fig = sprintf('%s/maaslin_ridges_%s_gender.tiff', out_path, taxa_sig_var)
# ggsave(out_fig,plot = p_ridge, device = 'tiff', width = 15,height = 10, compression  = 'lzw', bg = "white")
# 
# # maaslin2_plot_colored_scatters(taxa, x_var = taxa_sig_var, col_var = group_var, head_path, name, out_dir_name = sprintf('my_%s_scatters', taxa_sig_var), fdr = 0.25)
# 
# # p_heatmap25 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', taxa_sig_var = taxa_sig_var, fdr = 0.25, taxonomy_flag = F)
# # p_heatmap1 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', taxa_sig_var = taxa_sig_var, fdr = 0.1)
# # # taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)
# # # p_heatmap25 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_all_samples', taxa_sig_var = taxa_sig_var, fdr = 0.25)
# # # p_heatmap1 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_all_samples', taxa_sig_var = taxa_sig_var, fdr = 0.1)
# # # 
# 
# # path = '../res_israel_noRS_amnonFlt/'
# # taxa_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, 'source')
# # taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)
# # # p_heatmap25 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_israel_noRS', taxa_sig_var = taxa_sig_var, fdr = 0.25)
# # # p_heatmap1 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_israel_noRS', taxa_sig_var = taxa_sig_var, fdr = 0.1)
# 
# # fdr = 0.25
# # taxa_m = maaslin2_set_data_for_heatmap_plot(taxa, head_path, name, x_var, taxa_sig_var, group_var = group_var, taxonomy_file = taxonomy_file, fdr)
# # taxa_m = merge(taxa_m, taxa[,c('SampleID','location')], by = 'SampleID', all.x=T, sort = F)
# # 
# # library(viridis)
# # p_heatmap = ggplot(taxa_m, aes_string(x=x_var, y='taxa', fill = 'log10_RA')) + geom_tile() + 
# #   facet_grid(~Patient_group2+location, scales = 'free') + 
# #   theme_classic() + 
# #   ylab('ASV') + 
# #   scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
# #   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
# #   scale_fill_viridis() 
# # out_path = sprintf('%s/%s/res/%s',head_path,name, 'my_heatmaps')
# # dir.create(out_path)
# # out_fig = sprintf('%s/%s_sigTaxa_%s_%s_heatmap%s_fdr%s.tiff', out_path, taxa_sig_var, x_var, group_var, '_israel_noRS', fdr)
# # ggsave(out_fig,plot = p_heatmap, device = 'tiff', width = 12,height = 2+length(unique(taxa_m$ASV))/8, compression  = 'lzw')
# 
