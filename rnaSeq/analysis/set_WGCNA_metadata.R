# source('/pita/users/tzipi/code/R_figs/WGCNA_funcs.R')
add_omic = 'FFQ'

core_multimap = '../multiomics/data/CORE_screening_multiomics_map_v6.txt'
map_file = core_multimap
metadata_df = read.table(map_file, header = T, sep = '\t', na.strings = c('NA','','#N/A'), comment.char = '',quote = '')
metadata_df = metadata_df[!is.na(metadata_df$ID_rnaSeq_TI),]
row.names(metadata_df) = make.names(metadata_df$ID_rnaSeq_TI)


# add 16S metadata
add_file = '/pita/users/tzipi/projects/multiomics/CORE/16S/res_stool_TI_v3_amnonClean/biom/CORE_deblur_map_L7.2.txt'
add_df = read.table(add_file, header = T, sep = '\t', quote = '', comment.char = '')

for ( i in 1:dim(add_df)[1] ) # in SOURCE cases where there 2 rnaSeq runs, match samples
{
  if ( !is.na(add_df$ID_rnaSeq_TI[i]) & add_df$Cohort[i] == 'SOURCE' )
  {
    if ( !add_df$ID_rnaSeq_TI[i] %in% row.names(metadata_df) )
    {
      temp = row.names(metadata_df)[ metadata_df$pn_ID == add_df$pn_ID[i] ]
      if ( length(temp) == 1 )
        add_df$ID_rnaSeq_TI[i] = temp
    }
  }
}

taxa_vars = c('health_index','dysbiosis_index','faith_pd',
              'IBD_index','dbbact_saliva','dbbact_crohns','dbbact_rural','dbbact_monkey')
add_df = add_df[,c('ID_rnaSeq_TI','Source',taxa_vars)]
add_df_st = add_df[add_df$Source == 'stool',]
add_df_ti = add_df[add_df$Source == 'TI',]

names(add_df_st) = sprintf('%s_stool', names(add_df_st))
names(add_df_ti) = sprintf('%s_TI', names(add_df_ti))

metadata_df$ID_rnaSeq_TI = row.names(metadata_df)
metadata_df2 = merge(metadata_df, add_df_st, by.x = 'ID_rnaSeq_TI', by.y = 'ID_rnaSeq_TI_stool', all.x = T)
metadata_df2 = merge(metadata_df2, add_df_ti, by.x = 'ID_rnaSeq_TI', by.y = 'ID_rnaSeq_TI_TI', all.x = T)
metadata_df = metadata_df2
row.names(metadata_df) = metadata_df$ID_rnaSeq_TI

## add FFQ data
ffq_merge_file = core_multimap
# ffq_ids = read.table(ffq_merge_file, header = T, sep = '\t', na.strings = c('NA','','#N/A'), row.names = 1, comment.char = '',quote = '')
ffq_ids = read.table(ffq_merge_file, header = T, sep = '\t', na.strings = c('NA','','#N/A'), comment.char = '',quote = '')
# ffq_file = '/pita/users/tzipi/projects/multiomics/IBD_cohorts/FFQ/data/FFQ_dec24/FFQ_dec24_org_norm_corf.txt'
ffq_file = '../FFQ/data/CORE_FFQ.txt'
ffq = read.table(file = ffq_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'), quote = '', comment.char = '')
ffq = ffq[,!names(ffq) =='Tot_fat_per_kcal']
# maaslin_res_file = '../FFQ/res/maaslin2/Dx_activity_age_gender_ctlRef/res/significant_results.tsv'
# maas = read.table(file = maaslin_res_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'), quote = '', comment.char = '')
# ffq = ffq[,c('pn_ID',unique(maas$feature[maas$metadata == 'Dx_activity']))]
# add_name = sprintf('%s_maasDx', add_name)


add_df = ffq[,!names(ffq) %in% c('ID_FFQ','ffq_origin','Cohort','Dx',#'pn_ID',
                                 'Dx_cohort','first_sample')]
add_df = merge(x=ffq_ids[,c('pn_ID','ID_rnaSeq_TI')], y= add_df, by = 'pn_ID', all.y = T)
add_df = add_df[!is.na(add_df$ID_rnaSeq_TI),]

ffq_vars = names(add_df)[!names(add_df) %in% c('ID_FFQ','ID_rnaSeq_TI','pn_ID')]

metadata_df$ID_rnaSeq_TI = row.names(metadata_df)
metadata_df$ID_rnaSeq_TI = make.names(metadata_df$ID_rnaSeq_TI)
add_df$ID_rnaSeq_TI = make.names(add_df$ID_rnaSeq_TI)
metadata_df2 = merge(metadata_df, add_df, by = 'ID_rnaSeq_TI', all.x = T)
metadata_df = metadata_df2
row.names(metadata_df) = metadata_df$ID_rnaSeq_TI

# add_name = sprintf('%s_%s',add_name, 'FFQ')

## add metabolomics data
met_merge_file = core_multimap
# ffq_ids = read.table(ffq_merge_file, header = T, sep = '\t', na.strings = c('NA','','#N/A'), row.names = 1, comment.char = '',quote = '')
met_ids = read.table(met_merge_file, header = T, sep = '\t', na.strings = c('NA','','#N/A'), comment.char = '',quote = '')
# ffq_file = '/pita/users/tzipi/projects/multiomics/IBD_cohorts/FFQ/data/FFQ_dec24/FFQ_dec24_org_norm_corf.txt'
met_files = '../metabolomics/data/feces-mbx12-s106-f417_biom.txt'
met = read.table(file = met_files, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'), quote = '', comment.char = '', skip = 1, row.names = 1)
met = as.data.frame(t(met))

## normalize metabolomics data
added_data_type = 'metabolomicsLogClean'
mt = met
if( grepl('Log',added_data_type) & grepl('metabolomics',added_data_type) )
{
  if( grepl('LogClean',added_data_type) )
  {
    for ( i in 1:dim(mt)[2] )
    {
      # change zeroes to 5th of the smallest value per metabolite
      mt[mt[,i] == 0,i] = min(mt[mt[,i] != 0,i], na.rm = T)/5
      # change values over 4 SDs from the mean to the top/bottom value without it
      top_cut = mean(mt[,i], na.rm = T) + 4*sd(mt[,i], na.rm = T)
      mt[mt[,i] > top_cut,i ] = max(mt[ mt[,i]<=top_cut,i ], na.rm = T)
      bottom_cut = mean(mt[,i], na.rm = T) - 4*sd(mt[,i], na.rm = T)
      mt[mt[,i] < bottom_cut,i ] = min(mt[ mt[,i]>=bottom_cut,i ], na.rm = T)
    }
    # log transform
    mt = log10(mt)
  } else
    mt = log10(mt + 0.00001)
}
met= mt


met$ID_MBX_stool = row.names(met)
# maaslin_res_file = '../FFQ/res/maaslin2/Dx_activity_age_gender_ctlRef/res/significant_results.tsv'
# maas = read.table(file = maaslin_res_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'), quote = '', comment.char = '')
# ffq = ffq[,c('pn_ID',unique(maas$feature[maas$metadata == 'Dx_activity']))]
# add_name = sprintf('%s_maasDx', add_name)

add_df = met
met_ids$ID_MBX_stool = make.names( met_ids$ID_MBX_stool )
met_ids$ID_MBX_stool[met_ids$ID_MBX_stool == 'NA.']= NA
add_df = merge(x=met_ids[,c('ID_rnaSeq_TI','ID_MBX_stool')], y= add_df, by = 'ID_MBX_stool', all.x = T)
add_df = add_df[!is.na(add_df$ID_rnaSeq_TI),]

met_vars = names(add_df)[!names(add_df) %in% c('ID_MBX_stool','ID_rnaSeq_TI','pn_ID')]

metadata_df$ID_rnaSeq_TI = row.names(metadata_df)
metadata_df$ID_rnaSeq_TI = make.names(metadata_df$ID_rnaSeq_TI)
add_df$ID_rnaSeq_TI = make.names(add_df$ID_rnaSeq_TI)
metadata_df2 = merge(metadata_df, add_df, by = 'ID_rnaSeq_TI', all.x = T)
metadata_df = metadata_df2
row.names(metadata_df) = metadata_df$ID_rnaSeq_TI


metadata_df_old = metadata_df

if( add_omic == 'fecal_mets')
  add_vars = met_vars
if( add_omic == 'FFQ')
  add_vars = ffq_vars


wanted_vars = c('Dx',
                'Age','Gender','BMI',
                'Smoking_YN',
                # 'city_peripherality_index_2020_val',
                'city_socioeconomi_rank_2019',
                'CRP','FCP',
                # 'gastro_total_lewis_score',
                # 'gastro_highest_lewis_score',
                met_vars,
                ffq_vars,
                # 'health_index_stool',# 'health_index_TI',
                'dysbiosis_index_stool',# 'dysbiosis_index_TI',
                'faith_pd_stool',# 'faith_pd_TI',
                # 'IBD_index_stool','IBD_index_TI',
                # 'dbbact_saliva_stool',# 'dbbact_saliva_TI',
                # 'dbbact_crohns_stool','dbbact_crohns_TI',
                # 'dbbact_rural_stool','dbbact_rural_TI',
                # 'dbbact_monkey_stool','dbbact_monkey_TI',
                'Dx_activity','new_diagnosis')
                # 'IBD_Related_Surgeries_YN',
                # 'gastro_Previous_intestinal_surgery',
                # 'num_of_drugs',
                # 'any_bioligics_CD','any_bioligics_remission',
                # 'Adalimumab_remission',
                # 'No_treatment_or_5asa_CD','No_treatment_or_5asa_remission')
# 'gastro_Is_patient_currently_taking_any_of_the_following_medications._.check_all_that_apply._.choice.none.',
# 'gastro_Is_patient_currently_taking_any_of_the_following_medications._.check_all_that_apply._.choice.Imuran.6MP_.purinethol..',
# 'gastro_Is_patient_currently_taking_any_of_the_following_medications._.check_all_that_apply._.choice.Adalimumab_.Humira..',
# 'gastro_Is_patient_currently_taking_any_of_the_following_medications._.check_all_that_apply._.choice.Infliximab_.Remicade..',
# 'gastro_Is_patient_currently_taking_any_of_the_following_medications._.check_all_that_apply._.choice.Vedolizumab_.Entyvio..',
# 'gastro_Is_patient_currently_taking_any_of_the_following_medications._.check_all_that_apply._.choice.Ustekinumab_.Stelara..',
# 'gastro_IUS__Blood_flow_.Limbery_Score.',
# 'gastro_IUS__Maximal_SB_thickness_.mm.',
# 'gastro_IUS__Maximal_colon_thickness_.mm.',
# 'gastro_dass21_total_depression',
# 'gastro_dass21_total_anxiety',
# 'gastro_dass21_total_stress')
metadata_df = metadata_df[,wanted_vars]

# # set metadata
# metadata_df = metadata_df[,! names(metadata_df) %in% 
#                            c('Cohort','Dx_specific','Source','pn_ID')]
# #                               'Country','Library_Prep','ID2','Type_R',
# #                               'Sum_Study_Bx_Dx','Type_R_Bx','MicMac',
# #                               'SampleID_16S')]
metadata_df$Dx = ifelse(metadata_df$Dx == 'Control',0,1)
metadata_df$Gender = ifelse(metadata_df$Gender == 'Male',1,0)
metadata_df$Smoking_YN = ifelse(metadata_df$Smoking_YN == 'Yes',1,0)
# metadata_df$Active_vs_control = metadata_df$Dx
# metadata_df$Active_vs_control[metadata_df$Dx_activity == 'CD_remission'] = NA
metadata_df$Active_vs_remission = metadata_df$Dx
metadata_df$Active_vs_remission[metadata_df$Dx_activity == 'Control'] = NA
metadata_df$Active_vs_remission[metadata_df$Dx_activity == 'CD_remission'] = 0
metadata_df$Remission_vs_control = metadata_df$Dx
metadata_df$Remission_vs_control[metadata_df$Dx_activity == 'CD_Active'] = NA

# metadata_df$CD_new_vs_followup = metadata_df$Dx
# metadata_df$CD_new_vs_followup[metadata_df$new_diagnosis == 'FollowUp'] = 0
# metadata_df$CD_new_vs_followup[metadata_df$new_diagnosis == 'Control'] = NA
# metadata_df$CD_new_vs_control = metadata_df$Dx
# metadata_df$CD_new_vs_control[metadata_df$new_diagnosis == 'FollowUp'] = NA
# metadata_df$CD_followup_vs_control = metadata_df$Dx
# metadata_df$CD_followup_vs_control[metadata_df$new_diagnosis == 'Diagnosis'] = NA

# metadata_df$gastro_Previous_intestinal_surgery = ifelse(metadata_df$gastro_Previous_intestinal_surgery == 'Yes',1,0)

metadata_df$FCP[metadata_df$Dx == 0] = NA
metadata_df$CRP[metadata_df$Dx == 0] = NA

metadata_df = metadata_df[,! names(metadata_df) %in% 
                            c('Dx_activity','new_diagnosis')]

cat_df = data.frame(var = names(metadata_df), cat = NA, n = NA, group = NA)
cat_df$cat = 'General'
taxa_vars = sprintf('%s_stool',taxa_vars)
cat_df$cat[cat_df$var %in% taxa_vars] = 'Bacteria'
cat_df$cat[cat_df$var %in% ffq_vars] = 'FFQ'
cat_df$cat[cat_df$var %in% met_vars] = 'Metabolites'

for (i in 1:dim(cat_df)[1])
{
  pos = !is.na(metadata_df[[cat_df$var[i]]])
  cat_df$n[i] = sum( pos )
  temp = table(metadata_df_old$Dx_activity[ pos ] )
  cat_df$group[i] = paste(names(temp), temp, sep = ": ", collapse = ", ")
}

# metadata_df = metadata_df[, !names(metadata_df) == 'Dx']

write.table(metadata_df, 'data/WGCNA_map.txt', row.names = T, sep = '\t', quote = F)
write.table(cat_df, 'data/WGCNA_map_cats.txt', row.names = F, sep = '\t', quote = F)




