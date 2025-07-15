source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

norm_flag = T

na_str = c('no_data','_','NA','unknown', 'other','na','No_followup','ND','Unkown','',' ','.')

ffq_file = '../data/FFQ_dec24/FFQ_dec24_org_norm.txt'
ffq = read.table(file = ffq_file, header = T, sep = '\t', comment.char = '', quote = '', na.strings = na_str)


cor_ftr_file = '../data/FFQ_dec24/help_tables/correlated_FFQ_features_groups.txt'
cor_df = read.table(file = cor_ftr_file, header = T, sep = '\t', comment.char = '', quote = '', na.strings = na_str)

ffq_f = ffq[, !names(ffq) %in% cor_df$FFQ_feature[cor_df$what_to_do != 'as is']]

z_groups = unique(cor_df$corr_group[cor_df$what_to_do == 'z score average'])
for (zg in z_groups)
{
  fdf = ffq[,cor_df$FFQ_feature[cor_df$corr_group == zg]]
  fdf = as.data.frame( scale(fdf, center = T, scale = T) )
  zscor_mean = rowMeans(fdf)
  vars_names = gsub('_mcg.d_cal','', names(fdf))
  vars_names = gsub('_g.d_cal','', vars_names)
  vars_names = gsub('_serv.wk_cal','', vars_names)
  
  name = sprintf('%s_zmean', paste( vars_names, collapse = '+'))
  ffq_f[[name]] = zscor_mean
}
ffq_f = ffq_f[,!names(ffq_f) %in% c('Number_of_food_items','ffq_origin')]
ffq_f$Number_of_food_items = ffq$Number_of_food_items
ffq_f$ffq_origin = ffq$ffq_origin

## adding back fat percentage, as I need it
ffq_f$Tot_fat_per_kcal = ffq$Tot_fat_per_kcal

write.table(ffq_f, '../data/FFQ_dec24/FFQ_dec24_org_norm_corf.txt', 
            quote = F, sep = '\t', row.names = F)

cor_names_file = '../data/FFQ_dec24/help_tables/correlated_FFQ_features_names.txt'
ndf = read.table(file = cor_names_file, header = T, sep = '\t', comment.char = '', quote = '', na.strings = na_str)

names(ffq_f) = make.names(names(ffq_f))
for ( i in 1:dim(ndf)[1] )
{
  names(ffq_f)[names(ffq_f) == ndf$org_name[i]] = ndf$new_name[i]
}
names(ffq_f) = gsub('_cal$','_per_kcal',names(ffq_f) )
names(ffq_f) = gsub('Vegetables_.non.starchy.','Vegetables',names(ffq_f), fixed = T )
names(ffq_f) = gsub('Other_oil','Oils_excl_olive',names(ffq_f), fixed = T ) ## need to confirm

write.table(ffq_f, '../data/FFQ_dec24/FFQ_dec24_org_norm_corf_v2.txt', 
            quote = F, sep = '\t', row.names = F)


