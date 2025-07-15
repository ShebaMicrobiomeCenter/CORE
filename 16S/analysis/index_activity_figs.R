library(ggplot2)

source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

name = 'CORE'
if ( name == 'CORE')
{
  path = '../res_stool_TI_v3_amnonClean/'
  taxa_file = sprintf('%s/biom/%s_deblur_map_L7.2.txt',path, name)
  SV_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, name)
  a_file = sprintf('%s/core-metrics-results/faith_pd_vector/alpha-diversity.tsv',path)
  b_file = sprintf('%s/biom/core-metrics-results/unweighted_unifrac_distance_matrix/distance-matrix.tsv',path)
  pcoa_file = sprintf('%s/core-metrics-results/unweighted_unifrac_pcoa_results/ordination.txt',path)
  out_path =  sprintf('%s/R_res/',path)
}
na_str = c('no_data','_','NA','unknown', 'other','na','No_followup','')
taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str, comment.char = '', quote = '')
SV = read.table(SV_file,sep="\t", header=TRUE, na.strings = na_str, row.names = 1, comment.char = '', quote = '')

# out_path = sprintf('%s/activity_index_boxplots_all/',out_path)

# taxa = taxa[taxa$Source == 'TI',]
# out_path = sprintf('%s/basic_sigs_TI/',out_path)

taxa = taxa[taxa$Source == 'stool',]
out_path = sprintf('%s/activity_index_boxplots_stool/',out_path)

# dir.create(out_path)

# y_var = 'X_health_index'; y_lb = 'Health\nindex'
y_var = 'faith_pd2'; y_lb = 'Faith\'s\nPD'
# y_var = 'dbbact_saliva'; y_lb = 'dbBact\nsaliva'
# y_var = 'dysbiosis_index'; y_lb = 'Dysbiosis\nindex'
# y_var = 'IBD_index'; y_lb = 'IBD\nspecific'

check_val = 'Dx_activity'

health_res = plot_variable_and_stats_v2(taxa, check_val, y_val =y_var, y_lb = y_lb, sig_ast_flag = T, print_pvals = T)

taxa$Dx_activity[taxa$Dx_activity == 'CD_Active'] = 'CD active'
taxa$Dx_activity[taxa$Dx_activity == 'CD_remission'] = 'CD remission'

check_val='Dx_activity'
temp = taxa[[check_val]]
for ( i in 1:length(taxa[[check_val]] ) )
{
  taxa[[check_val]][i] = sprintf('%s n=%d',temp[i], sum( temp == temp[i], na.rm = T ) )
  if (is.na(temp[i]))
    taxa[[check_val]][i] = sprintf('%s n=%d',temp[i],  sum( is.na(temp) ) )
}

cols = c('#ff8080','#ffc080','#8d8dff')
g= ggplot(taxa) + 
  geom_boxplot(aes_string(x=check_val, y=y_var, fill = check_val)) + 
  theme_bw() + 
  scale_fill_manual(values = cols ) +
  scale_x_discrete( breaks = sort(unique(taxa[[check_val]])),
                    labels = gsub('.* ','', sort(unique(taxa[[check_val]]))) ) +
  theme(panel.grid =element_blank(), 
        legend.position = 'none',
        axis.text = element_text(size=10, colour = 'black'),
        axis.text.x = element_text(size=6, colour = 'black'),
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  xlab(y_lb) +ylab('')

res = add_significant_asterix_to_plot_BH_v2(p = g, var = as.factor(taxa[[check_val]]), 
                                            val = taxa[[y_var]], test_type = 'wilcox', 
                                            angle = 90, asterix_scale_var = 2.1, label_size = 4, line_width = 0.5)
g=res[[1]]
g = g+coord_flip()
out_file = sprintf( '%s/%s.tiff', out_path, sprintf('%s_%s',y_var, check_val) )
# ggsave(out_file,plot = g, device = 'tiff', width = 4,height = 1, compression = 'lzw')
ggsave(out_file,plot = g, device = 'tiff', width =3.5,height = 1, compression = 'lzw')

