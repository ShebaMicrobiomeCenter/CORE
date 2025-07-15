source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

norm_flag = T

na_str = c('no_data','_','NA','unknown', 'other','na','No_followup','ND','Unkown','',' ','.')

ffq_file = '../data/FFQ_dec24/iace_ExSurv_25122024_xls_FFQ_from_kathy.txt'
ffq = read.table(file = ffq_file, header = T, sep = '\t', comment.char = '', quote = '', na.strings = na_str)

ffq_old_file = '../data/old/iace_ExSurv_16nov24_excel_only.txt'
ffq_old = read.table(file = ffq_old_file, header = T, sep = '\t', comment.char = '', quote = '', na.strings = na_str)


ffq = ffq[,!names(ffq) %in% c('MUFA_cis_g.d','PUFA_cis_g.d')]
names(ffq) = gsub('__','',names(ffq))

## remove problematic samples: 
# SK0729 manually fixed in the sheets already, will be included next time. 
# CORE_074_3M_2 is almost the same as CORE_074_3M_1
ffq = ffq[!ffq$pn_ID %in% c('SK0729','CORE_074_3M_2'),]

## remove low calories etc.
ffq = ffq[ffq$Energy_kcal.d >= 800 & !is.na(ffq$Energy_kcal.d), ]
ffq = ffq[ffq$Water_g.d >= 500, ]


# # remove the medias components (such as Olive_oil_more_than_other_oils)
medas_cols = c()
for ( i in 1:dim(ffq)[2] )
  if ( all(ffq[,i] %in% c(0,1,NA)) )
    medas_cols = c(medas_cols, names(ffq)[i] )
ffq = ffq[, !names(ffq) %in% medas_cols]

if (norm_flag)
{
  ## norm by calories
  first_ffq = 'Water_g.d'
  last_ffq = 'Fiber_g.1000kcal'
  res = norm_ffq_by_energy_v2(ffq, first_ffq = first_ffq,
                           last_ffq = last_ffq,
                           ffq_norm_var = 'Energy_kcal.d',
                           no_norm_var_part = 'kcal',
                           no_norm_vars = c('Number_of_food_items',
                                            'NOVA_class','I.MEDAS_SCORE'),
                           norm_suffix = '_cal')
  ffq = res[[1]]
}


food_prms = names(ffq)
food_prms = food_prms[!food_prms %in% c('pn_ID','Source_of_FFQ_data')]

ffq_f = ffq[,food_prms]
row.names(ffq_f) = ffq$pn_ID

# bad_pos = which(colSums(ffq_f, na.rm = T)==0)
# ffq_f = ffq_f[,-bad_pos]

## cehck origin differnces. no problem in dec24 version
bad_ffq = c()
for (i in 1:dim(ffq_f)[2])
{
  pos = ffq$Source_of_FFQ_data == '1st Gastro-SOURCE'
  temp = mean(ffq_f[pos,i], na.rm = T) / mean(ffq_f[!pos,i], na.rm = T)
  if (is.finite(temp))
  {
    if (temp < 1)
      temp = mean(ffq_f[! pos,i], na.rm = T) / mean(ffq_f[pos,i], na.rm = T)
    if ( temp > 5 ) 
    {
      # print(sprintf('%s FC %s', names(ffq_f)[i], temp))
      bad_ffq = c(bad_ffq, names(ffq_f)[i])
    }
  }
}
# ffq_f[ffq$origin == 'laptop', names(ffq_f) %in% bad_ffq] = NA
# ffq_f = ffq_f[,!names(ffq_f) %in% bad_ffq]

## handle outliers
sd_cutoff = 5
for ( i in 1:dim(ffq_f)[2] )
{
  ##  change zeroes to 5th of the smallest value per metabolite
  # ffq_f[met[,i] == 0,i] = min(ffq_f[ffq_f[,i] != 0,i], na.rm = T)/5
  # change values over X SDs from the mean to the top/bottom value without it
  top_cut = mean(ffq_f[,i], na.rm = T) + sd_cutoff*sd(ffq_f[,i], na.rm = T)
  pos = ffq_f[,i] > top_cut & !is.na(ffq_f[,i] > top_cut)
  if ( sum(pos) > 0)
  {
    print(names(ffq_f)[i])
    print(i)
    ffq_f[pos,i ] = max(ffq_f[ !pos,i ], na.rm = T)
  }
}




ffq_f2 = data.frame(pn_ID = ffq$pn_ID, ffq_f, ffq_origin = ffq$Source_of_FFQ_data)
ffq = ffq_f2

## add metdata
temp = data.frame(pn_ID = ffq$pn_ID)
pos = grepl('CORE',ffq$pn_ID)
ffq$ffq_origin[pos] = sprintf('%s_CORE', ffq$ffq_origin[pos])
temp$Cohort = ifelse(pos, 'CORE','seker')
temp$Cohort[grepl('SOURCE',ffq$pn_ID)] = 'SOURCE'
temp$Cohort[grepl('CURE',ffq$pn_ID)] = 'CURE'
temp$Cohort[grepl('IIRN',ffq$pn_ID)] = 'IIRN'
temp$Dx = 'Control'
temp$Dx[temp$Cohort %in% c('CORE','CURE','IIRN')] = 'CD'
source_file = '/pita/users/tzipi/projects/multiomics/SOURCE/metadata/SOURCE_multiomics_map.txt'
source_df = read.table(file = source_file, header = T, sep = '\t', comment.char = '', quote = '', na.strings = na_str)
source_df$pn_ID = sprintf('SOURCE_%s',source_df$pn_ID)
temp$Dx[ffq$pn_ID %in% source_df$pn_ID[source_df$Patient_group2 == 'Israel_CD']] = 'CD'
temp$Dx[ temp$Cohort == 'SOURCE' & ! ffq$pn_ID %in% source_df$pn_ID] = NA
# manually add missing based on the IBD full Dx table
temp$Dx[temp$pn_ID %in% c('SOURCE_A049')] = 'CD'
temp$Dx[temp$pn_ID %in% c('SOURCE_A048','SOURCE_A037')] = 'UC'

temp$Dx_cohort = sprintf('%s_%s',temp$Dx, temp$Cohort)

## seperate sampleID and patient ID
temp$ID = temp$pn_ID
temp$ID = gsub('_[1-9]$','',temp$ID)
temp$ID = gsub(' Copy$','',temp$ID)
temp$ID = gsub('_[0-9][0-9]_[0-9][0-9]_20[0-9][0-9]$','',temp$ID)
temp$ID = gsub(' .*$','',temp$ID)
temp$ID = gsub('_S$','',temp$ID)
temp$ID = gsub('_[0-9]+M$','',temp$ID)

## mark first sample in repeats. note this is manual
## check with updates!!
temp$first_sample = 'Yes'
reps = as.data.frame(table(temp$ID))
temp$first_sample[temp$ID %in% reps$Var1[reps$Freq > 1]] = 'No'
temp$first_sample[temp$first_sample == 'No' & grepl('_S$',temp$pn_ID)] = 'Yes'
temp$first_sample[temp$pn_ID %in% c('CORE_024_3M','CORE_022_3M')] = 'Yes'

ffq = merge(temp, ffq, by='pn_ID')
names(ffq)[names(ffq) == 'pn_ID'] = 'SampleID'
names(ffq)[names(ffq) == 'ID'] = 'pn_ID'


if (norm_flag)
{
  write.table(ffq, '../data/FFQ_dec24/FFQ_dec24_org_norm.txt', sep = '\t',row.names = F, quote = F)
} else
  write.table(ffq, '../data/FFQ_dec24/FFQ_dec24_org.txt', sep = '\t',row.names = F, quote = F)

