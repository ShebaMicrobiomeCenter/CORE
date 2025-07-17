library(ggplot2)

## rnaSeq
core_screening_file = '../data/CORE_screening_multiomics_map_v8.txt'
map = read.table(file = core_screening_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'), quote = '', comment.char = '')

# map$Medications = map$gastro_Is_patient_currently_taking_any_of_the_following_medications._.check_all_that_apply._.choice.none.
# map$Medications = 1-map$Medications
# map$Medications[map$Cohort %in% c('SOURCE','seker')] = 0
# map$Medications[map$new_diagnosis %in% c('Control','Diagnosis')] = 0
# map$Medications[is.na(map$Medications)] = 1 # for now, ask yael
# map$Medications[map$Medications == 0] = NA

map$Biologics = map$any_bioligics_CD
map$Biologics[map$any_bioligics_CD == 0] = NA

# map$active = map$Dx_activity
# map$active[map$Dx_activity != 'CD_Active'] = NA
map$New_diagnosis = map$new_diagnosis
map$New_diagnosis[map$New_diagnosis == 'CD_Active'] = 1
map$New_diagnosis[map$New_diagnosis == 'FollowUp'] = NA
map$New_diagnosis[map$New_diagnosis == 'Control'] = NA


map$highest_lewis_score_over_350 = 
  ifelse(map$gastro_highest_lewis_score > 350, 1, NA)

met_vars = c('pn_ID','Cohort','Dx','Dx_activity','new_diagnosis')
show_vars = c('New_diagnosis','Biologics','highest_lewis_score_over_350',# 'gastro_dass21_total_stress',
              'CRP','FCP','gastro_highest_lewis_score',
              'FFQ_available',
              'ID_16s_stool',# 'ID_16s_TI',
              'ID_MBX_stool','ID_MBX_stool_other_batch',
              # 'ID_MBX_serum','ID_MBX_serum_other_batch',
              'ID_rnaSeq_TI')
df = map[,c(met_vars, show_vars)]

df$FFQ_available[df$FFQ_available == 'no'] = NA
# df$CRP[df$Dx == 'Control'] = NA
df$FCP[df$Dx == 'Control'] = NA

df[,show_vars] = !is.na(df[,show_vars])

df$ID_MBX_stool[!df$ID_MBX_stool & df$ID_MBX_stool_other_batch  ] = 0.5
# df$ID_MBX_serum[!df$ID_MBX_serum & df$ID_MBX_serum_other_batch  ] = 0.5

df =df[,!names(df) %in% c('ID_MBX_stool_other_batch','ID_MBX_serum_other_batch')]
show_vars = show_vars[!show_vars %in% c('ID_MBX_stool_other_batch','ID_MBX_serum_other_batch')]
df$fill_sum = rowSums(df[,show_vars])

all(df$pn_ID == map$pn_ID)
# df$active[df$active == F] = 0.2
# df$active[df$Dx == 'Control'] = 0
df$New_diagnosis[df$New_diagnosis == F] = 0.2
df$New_diagnosis[df$Dx == 'Control'] = 0
# df$Medications[df$Medications == F] = 0.2
df$Biologics[df$Biologics == F] = 0.2
df$Biologics[is.na(map$any_bioligics_CD)] = 0

df$highest_lewis_score_over_350[df$highest_lewis_score_over_350 == F] = 0.2
df$highest_lewis_score_over_350[!df$gastro_highest_lewis_score] = 0

dfm = reshape::melt(data = df, id.vars = c(met_vars,'fill_sum'), variable_name = 'Data_type')

dfm$pn_ID = factor(dfm$pn_ID, levels = df$pn_ID[order(df$fill_sum)])

dfm$Patient_group3 = dfm$Dx_activity
for (pg in unique(dfm$Dx_activity))
{
  dfm$Patient_group3[dfm$Dx_activity == pg] = 
    sprintf('%s (n=%s)', gsub('_',' ',pg), 
            length(unique( dfm$pn_ID[dfm$Dx_activity == pg] )))
}

# dfm = dfm[dfm$Type !='metabolomics',]
dfm$Type = 'Metadata'
dfm$Type[dfm$Data_type %in% c('ID_16s_TI','ID_16s_stool')] = 'Bacteria'
dfm$Type[dfm$Data_type %in% c('ID_rnaSeq_TI')] = 'RNA-Seq'
dfm$Type[dfm$Data_type %in% c('FFQ_available')] = 'FFQ'
dfm$Type[dfm$Data_type %in% c('ID_MBX_serum','ID_MBX_stool')] = 'Metabolomics'
dfm$Type[dfm$Data_type %in% c('active','Medications','Biologics','New_diagnosis',
                              'highest_lewis_score_over_350')] = 'Disease characteristics'

dfm$Data_type = as.character(dfm$Data_type)
dfm$Data_type[dfm$Data_type == 'FFQ_available'] = 'FFQ'
dfm$Data_type[dfm$Data_type == 'ID_16s_stool'] = '16S Stool'
dfm$Data_type[dfm$Data_type == 'ID_16s_TI'] = '16S TI'
dfm$Data_type[dfm$Data_type == 'ID_MBX_stool'] = 'MBX Stool'
dfm$Data_type[dfm$Data_type == 'ID_MBX_serum'] = 'MBX Serum'
dfm$Data_type[dfm$Data_type == 'ID_rnaSeq_TI'] = 'RNA-Seq TI'
dfm$Data_type[dfm$Data_type == 'active'] = 'Disease activity Y/N'
dfm$Data_type[dfm$Data_type == 'New_diagnosis'] = 'Newly diagnosed Y/N'
dfm$Data_type[dfm$Data_type == 'Medications'] = 'Medications Y/N'
dfm$Data_type[dfm$Data_type == 'Biologics'] = 'Biologics Y/N'
dfm$Data_type[dfm$Data_type == 'gastro_highest_lewis_score'] = 'Lewis score'
dfm$Data_type[dfm$Data_type == 'highest_lewis_score_over_350'] = 'Lewis score > 350'
dfm$Data_type[dfm$Data_type == 'gastro_dass21_total_stress'] = 'DASS-21'


dfm$Data_type = factor(dfm$Data_type, levels = c(# 'Disease activity Y/N',
                                                 'Newly diagnosed Y/N',
                                                 'Biologics Y/N',
                                                 'Lewis score > 350',
                                                 # 'DASS-21',
                                                 'Lewis score','FCP','CRP','FFQ',
                                                 '16S Stool','16S TI',
                                                 'MBX Stool',# 'MBX Serum',
                                                 'RNA-Seq TI'))

dfm$Type = factor(dfm$Type, levels = c('Disease characteristics','Metadata','FFQ','Bacteria','Metabolomics','RNA-Seq'))
# dfm$Patient_group3 = factor(dfm$Patient_group3, 
#                             levels = c('FollowUp (n=85)','Diagnosis (n=31)',
#                                        'Control (n=77)'))
dfm$Patient_group3 = factor(dfm$Patient_group3, 
                            levels = c('CD Active (n=37)','CD remission (n=77)',
                                       'Control (n=77)'))

cols = c('gray20','#293462','#006E7F','#F8CB2E','#EE5007','#990000')
g = ggplot(dfm, aes(x=Data_type, y=pn_ID, alpha = as.character(value), fill = Type)) + 
  geom_tile() + 
  # geom_tile(colour = 'gray40') + 
  scale_alpha_manual(values = c(.9,0.6,0.4, 0), name = 'Data\navailable', 
                     breaks = c('1','0.5','0.2','0'), labels = c('Yes','Other batch','No','No')) +
  scale_fill_manual(values = cols) + 
  # facet_grid(~Type + Source, scales = 'free') + 
  coord_flip() + 
  facet_grid(.~Patient_group3 , scales = 'free', space = 'free') + 
  ylab('Subject') + xlab ('')  +
  theme_bw() +
  scale_y_discrete(expand=c(0,0)) +scale_x_discrete(expand=c(0,0)) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = "black"))

# # ggsave('Israel_cohort_figure_full.pdf',plot = g, device = 'pdf', width = 4,height = 8, limitsize = FALSE)
ggsave('CORE_cohort_figure_S_full.pdf',plot = g, device = 'pdf', width = 22,height = 4.5, limitsize = FALSE)
# 
g2 = g+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank() )
# ggsave('Israel_cohort_figure.pdf',plot = g2, device = 'pdf', width = 4,height = 8, limitsize = FALSE)
ggsave('CORE_cohort_figure_S.pdf',plot = g2, device = 'pdf', width = 10,height = 2.5, limitsize = FALSE)
write.table(x = df, file = 'CORE_cohort_figure_S.txt', sep = '\t', row.names = F, quote = F)

# 
# row.names(dx_map) = dx_map$pn_ID
# dx_map = dx_map[df$pn_ID,]
# df$Age = dx_map$Age_years
# df$Gender = dx_map$Gender
# df$BMI = dx_map$BMI
# df$CRP_mg_L = dx_map$CRP_numeric
# df$Calprotectin_ug_g = dx_map$Calprotectin_numeric
# 
# write.table(x = df, file = 'Israel_basic_map.txt', quote = F, sep = '\t',row.names = F)
# 
# 
