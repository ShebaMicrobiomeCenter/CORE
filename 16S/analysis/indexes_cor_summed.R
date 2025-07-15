library(ggplot2)
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')


name = 'CORE'
if ( name == 'CORE')
{
  path = '../res_stool_TI_v3_amnonClean/'
  # path = '../res_stool_v3_amnonClean/'
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

SV$faith_pd = SV$faith_pd2

ffq_file = '../../FFQ/data/CORE_FFQ_with_metadata.txt'
ffq = read.table(ffq_file, header = T, sep = '\t', row.names = 1)

multi_file = '../../multiomics/data/CORE_screening_multiomics_map_v8.txt'
multi_map = read.table(multi_file, header = T, sep = '\t', na.strings = c('NA','','#N/A'), comment.char = '',quote = '')
multi_map = multi_map[, c('ID_16s_stool', 'ID_FFQ','pn_ID')]
multi_map = multi_map[!is.na(multi_map$ID_16s_stool) & !is.na(multi_map$ID_FFQ),]

SV = SV[multi_map$ID_16s_stool,]
ffq = ffq[multi_map$pn_ID,]

SV = cbind(SV, ffq)

# categories_file = '../data/trasnlate_and_org_tables_nov24/seker_map_categories_feb25.txt'
# cat = read.table(categories_file, header = T, sep = '\t')
# cat$Category[cat$include == 0] = 'no'

# health_vars = c('health_index', 'faith_pd','IBD_index','dysbiosis_index',
#                 'dbbact_saliva','dbbact_crohns','dbbact_rural','dbbact_monkey')
# health_vars = c('health_index', 'faith_pd','IBD_index','dysbiosis_index')
health_vars = c( 'faith_pd','dysbiosis_index')


last_met = 'Description'
# last_met = 'remark'
# metadata_pos = 1:which(names(SV) == last_met)

metadata_pos = c( names(ffq), health_vars )
SV_map = SV[metadata_pos]

taxa_pos = ( which(names(SV) == last_met) +1 ):dim(SV)[2]
SV_RA = SV[,taxa_pos]

# SV_rank = RA_biom_to_rank(SV_RA)
# final_SV = SV_rank
# # final_SV = SV_RA
# # uniDI = calculate_health_index(final_SV, up_df, down_df)
# binary_health_index = calculate_binary_health_index(SV_RA, up_df, down_df)

# out_path_health = sprintf('%s/health_index_freq_sum/', out_path)
# dir.create(out_path_health)

taxa = SV_map


# metadata_cols = names(taxa)[1:which( names(taxa) == 'Description' )]
metadata_cols = names(taxa)
bad_vars = c('pn_ID', 'Age','BMI',
             health_vars)
metadata_cols = metadata_cols[! metadata_cols%in% bad_vars]

calc_cor_df = function(taxa, metadata_cols,health_var, type )
{
  res_df = data.frame(var = metadata_cols, 
                      spearman_r=NA, spearman_p=NA, n=NA, cor_var = health_var)
  for ( i in 1:length(metadata_cols) )
  {
    check_val = metadata_cols[i]
    ## check numeric and more than 50 values
    if ( is.numeric( taxa[[check_val]] ) & sum(!is.na(taxa[[check_val]])) > 15 )
    {
      ## if numeric that is converted from categoric (0,1 for example)
      ## make sure there are at least 2 categories with 25 samples each
      if ( ! ( length(unique(taxa[[check_val]] )) <= 5 & 
               sum(table(taxa[[check_val]]) >= 10 ) < 2 ) )
      {
        res = cor.test(taxa[[check_val]], taxa[[health_var]], method = 'spearman')
        res_df$spearman_r[i] = res$estimate
        res_df$spearman_p[i] = res$p.value
        res_df$n[i] = sum(!is.na(taxa[[check_val]]))
      }
    }
  }
  res_df$type = sprintf('%s (n=%s)', type, dim(taxa)[1])
  res_df = res_df[!is.na(res_df$spearman_r),]
  res_df$spearman_q = p.adjust(res_df$spearman_p, method = 'BH')
  
  return(res_df)
}

calc_cor_df_multiple = function(taxa, metadata_cols,health_vars, type )
{
  res_df = calc_cor_df(taxa, metadata_cols,health_vars[1], type )
  if (length(health_vars) > 1)
  {
    for (i in 2:length(health_vars))
    {
      temp = calc_cor_df(taxa, metadata_cols,health_vars[i], type )
      res_df = rbind(res_df, temp)
    }
  }
  return(res_df)
}


res_df = calc_cor_df_multiple(taxa, metadata_cols,health_vars, 'All' )

# temp = calc_cor_df_multiple(taxa[taxa$Gender == 'Male',], metadata_cols,health_vars, 'Male' )
# res_df = rbind(res_df, temp)
# 
# temp = calc_cor_df_multiple(taxa[taxa$Gender == 'Female',], metadata_cols,health_vars, 'Female' )
# res_df = rbind(res_df, temp)

temp = calc_cor_df_multiple(taxa[taxa$Dx_activity == 'CD_Active',], metadata_cols,health_vars, 'CD active' )
res_df = rbind(res_df, temp)

temp = calc_cor_df_multiple(taxa[taxa$Dx_activity == 'CD_remission',], metadata_cols,health_vars, 'CD remission' )
res_df = rbind(res_df, temp)

temp = calc_cor_df_multiple(taxa[taxa$Dx_activity == 'Control',], metadata_cols,health_vars, 'Control' )
res_df = rbind(res_df, temp)

temp = calc_cor_df_multiple(taxa[taxa$Dx == 'CD',], metadata_cols,health_vars, 'CD' )
res_df = rbind(res_df, temp)

temp = calc_cor_df_multiple(taxa[taxa$Dx == 'CD' & SV$FCP > 250,], metadata_cols,health_vars, 'CD FCP>250' )
res_df = rbind(res_df, temp)

temp = calc_cor_df_multiple(taxa[taxa$Dx == 'CD' & SV$FCP <= 250,], metadata_cols,health_vars, 'CD FCP<=250' )
res_df = rbind(res_df, temp)


res_df_temp = res_df
res_df_temp = res_df_temp[res_df_temp$type == 'All (n=119)',]
res_df_temp$spearman_r = abs(res_df_temp$spearman_r)
summary(res_df_temp$spearman_r[res_df_temp$cor_var == 'dysbiosis_index'])
summary(res_df_temp$spearman_r[res_df_temp$cor_var == 'faith_pd'])
table(res_df_temp$spearman_q[res_df_temp$cor_var == 'dysbiosis_index'] <= 0.25)
table(res_df_temp$spearman_q[res_df_temp$cor_var == 'faith_pd'] <= 0.25)

good_vars = unique( res_df$var[res_df$spearman_q < 0.25] )
# good_vars = unique( res_df$var[res_df$spearman_p < 0.05] )

res_df_all = res_df
res_df = res_df[res_df$var %in% good_vars, ]

res_df$type = factor(res_df$type, levels = c('All (n=119)',
                                             'CD (n=63)',
                                             'CD active (n=16)',
                                             'CD remission (n=47)',
                                             'CD FCP>250 (n=34)',
                                             'CD FCP<=250 (n=29)',
                                             'Control (n=56)'))
res_dft = res_df
res_dft$spearman_r[res_dft$cor_var == 'faith_pd'] = -1*res_dft$spearman_r[res_dft$cor_var == 'faith_pd']
var_mean_cor = aggregate(res_dft$spearman_r, list(res_dft$var), FUN=mean)
res_df$var = factor(res_df$var, var_mean_cor$Group.1[order(var_mean_cor$x)])

res_df$cor_var = factor(res_df$cor_var, levels = health_vars)


max_col = max(abs(res_df$spearman_r))

# res_df= merge(res_df, cat, by.x = 'var', by.y = 'Var', all.x = T)

col2 = colorRampPalette(rev( c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                               "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") ))(200)
res_df$cor_var = as.character(res_df$cor_var)
res_df$cor_var[res_df$cor_var == 'faith_pd'] = 'Faith\'s\nPD'
res_df$cor_var[res_df$cor_var == 'dysbiosis_index'] = 'Dysbiosis\nIndex'
res_df$cor_var = factor(res_df$cor_var, levels = c('Dysbiosis\nIndex','Faith\'s\nPD'))

cor_fig = ggplot(res_df, 
                 aes(x=type, y=var, colour = spearman_r, 
                     size= abs(spearman_r), 
                     label = ifelse(spearman_q<=0.1, '*','') )) + 
  geom_point(alpha = 1) + 
  geom_point(aes(shape = spearman_q<=0.25), size=1.5, colour = 'gray10', alpha=1) + 
  # geom_point(aes(shape = spearman_p<=0.05), size=1.5, colour = 'black', alpha= 0.6) +
  scale_shape_manual(values = c(NA, 3)) +
  # scale_fill_gradientn(colours = col2) +
  facet_grid(.~cor_var,scales = 'free', space = 'free', switch="y") + 
  scale_colour_gradientn(# colors=c("red","yellow","green"), 
    colors = col2,
    limits = c(-1*max_col,max_col), 
    name = 'Spearman\'s\nrho', na.value = 'gray') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        # panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        panel.grid = element_blank(),
        strip.text.y.left = element_text(angle = 0)) + 
  # geom_text() + 
  # geom_text(aes(label=sprintf('%.2f\n%.2e',value,q_val )), size=3)
  # scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  # geom_text(aes(label=sprintf('%.2f',value )), size=3) + 
  xlab('') + ylab('')
cor_fig

ggsave(sprintf('%s/cor_light_sum_q.pdf', out_path),
       plot = cor_fig,  device = 'pdf',
       # width = 6.2,height = 8.5)
       width = 6.4,height = 8.5)
       # width = 10,height = 10)
       # width = 16,height = 10)
       # width = 13,height = 20)


