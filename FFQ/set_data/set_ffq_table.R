source('/pita/users/tzipi/code/R_figs/figs_funcs.R')


## filter full FFQ file to samples relevant to the CORE cohort
core_file = '../../multiomics/data/CORE_screening_multiomics_map_v8.txt'
# ffq_ids = read.table(ffq_merge_file, header = T, sep = '\t', na.strings = c('NA','','#N/A'), row.names = 1, comment.char = '',quote = '')
core_ids = read.table(core_file, header = T, sep = '\t', na.strings = c('NA','','#N/A'), comment.char = '',quote = '')
# ffq_file = '/pita/users/tzipi/projects/multiomics/IBD_cohorts/FFQ/data/FFQ_dec24/FFQ_dec24_org_norm.txt'
# ffq_file = '/pita/users/tzipi/projects/multiomics/IBD_cohorts/FFQ/data/FFQ_dec24/FFQ_dec24_org_norm_corf.txt'
ffq_file = '/pita/users/tzipi/projects/multiomics/IBD_cohorts/FFQ/data/FFQ_dec24/FFQ_dec24_org_norm_corf_v2.txt'
ffq = read.table(file = ffq_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'), quote = '', comment.char = '')
ffq$pn_ID = gsub('CORE_','CORE',ffq$pn_ID)

ffq = ffq[,!names(ffq) == 'Maganese_mg.d_cal']

core_ids = core_ids[!is.na(core_ids$ID_FFQ),]

dff = ffq[,!names(ffq) %in% c('ffq_origin','Cohort','Dx','Dx_cohort','first_sample')]
rownames(dff) = dff$SampleID
dff = dff[core_ids$ID_FFQ, ]
# use same names for pn_ID
# because we used some same patient in different cohort samples
# also to match for metadata in the future
dff$pn_ID = core_ids$pn_ID 

dff = dff[dff$SampleID %in% core_ids$ID_FFQ, ]

dff = dff[,!names(dff) == 'SampleID']


## handle outliers
sd_cutoff = 5
for ( i in 1:dim(dff)[2] )
{
  ##  change zeroes to 5th of the smallest value per metabolite
  # ffq_f[met[,i] == 0,i] = min(ffq_f[ffq_f[,i] != 0,i], na.rm = T)/5
  # change values over X SDs from the mean to the top/bottom value without it
  top_cut = mean(dff[,i], na.rm = T) + sd_cutoff*sd(dff[,i], na.rm = T)
  pos = dff[,i] > top_cut & !is.na(dff[,i] > top_cut)
  if ( sum(pos) > 0)
  {
    print(names(dff)[i])
    print(i)
    dff[pos,i ] = max(dff[ !pos,i ], na.rm = T)
  }
}


## fix weird values
## REMOVE IN NEWER VERISONS, THIS IS SUPER SPECIFIC!!
# dff$Maganese_mg.d_cal[dff$Maganese_mg.d_cal > 20] = 7
# dff$Biotin_mcg.d_cal[dff$Biotin_mcg.d_cal > 18] = 4

write.table(dff, '../data/CORE_FFQ.txt', quote = F, sep = '\t', row.names = F)

mets_vars = c('Cohort','Dx','Dx_activity','new_diagnosis',
              'Age','Gender','BMI','Smoking_YN')

dff2 = merge(x=core_ids[c('pn_ID',mets_vars)], y=dff, by = 'pn_ID' , all.y = T)
write.table(dff2, '../data/CORE_FFQ_with_metadata.txt', quote = F, sep = '\t', row.names = F)

dff3 = merge(x=core_ids, y=dff, by = 'pn_ID' , all.y = T)
write.table(dff3, '../data/CORE_FFQ_with_all_metadata.txt', quote = F, sep = '\t', row.names = F)

