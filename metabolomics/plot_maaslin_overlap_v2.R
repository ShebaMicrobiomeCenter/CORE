source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(ggplot2)

na_str = c('no_data','_','NA','unknown', 'other','na','No_followup')

core_maas_file = 'maasin2_res/CORE_Feces_Dx_gender_age/res/significant_results.tsv'
mc = read.table(core_maas_file,sep="\t", header=TRUE, na.strings = na_str, quote = '',comment.char = '')

source_maas_file = 'SOURCE_maaslin_all_results.tsv'
ms = read.table(source_maas_file,sep="\t", header=TRUE, na.strings = na_str, quote = '',comment.char = '')

mc = mc[mc$metadata == 'Dx',]
ms = ms[ms$metadata == 'Dx',]
ms = ms[ms$feature %in% mc$feature,]
row.names(mc) = mc$feature
mc = mc[ms$feature,]
names_mc = c('feature','metadata','value','coef','qval')
new_names_mc = c('feature','metadata','value','coef_core','qval_core')
mc = mc[,names_mc]
names(mc) = new_names_mc
mc$coeff_source = ms$coef
mc$q.value_source = ms$qval

df = mc
# file_path = 'data/maas_overlap_yael.txt'
# df = read.table(file_path,sep="\t", header=TRUE, na.strings = na_str, quote = '',comment.char = '')

df$coef_core = -1*df$coef_core 
df$coeff_source = -1*df$coeff_source 

# dif_text = 'Not significant in\nCD active vs control'
# dif_text = 'Not significant'
dif_text = 'Different directions'

Direction = ifelse(df$coef_core < 0,'Higher in controls','Higher in CD')
Direction[ (df$coef_core > 0) != (df$coeff_source > 0) ] = dif_text
Feature2 = df$feature
# remove long names
Feature2[Feature2 %in% c(# 'Ketoleucine',
  'all.cis.4.7.10.13.16.Docosapentaenoic.acid')] = NA
Feature2[Direction == dif_text] = NA
library(ggrepel)
g = ggplot( df, 
            # aes(x=coef_rural, y=coef_dx, colour = q_dx<=0.25) ) + 
            aes(x=coef_core, y=coeff_source, label = Feature2) ) +
  geom_vline(xintercept = 0, colour ='gray80', linetype="dashed") + 
  geom_hline(yintercept = 0, colour ='gray80', linetype="dashed") + 
  scale_colour_manual(values = c('gray60','#c10000','#0070c0'), name = '') +
  # scale_colour_manual(values = c('#c10000','#0070c0','gray60'), name = '') +
  # geom_text_repel(size = 2, seed = 4, nudge_y = 0.2, box.padding = 7) + 
  geom_point(aes(colour = Direction)) + theme_bw() + theme(panel.grid = element_blank()) + 
  geom_text_repel(size = 2, seed=2023,  box.padding = 7) +
  # aes(x=FC_rural, y=FC_dx, colour = q_dx<=0.25) ) + 
  xlab('CORE\nCD remission vs control') + ylab('SOURCE\nCD active vs control') 
# ggsave(sprintf('res/source_core_maas_overlap_scatter.tiff'),
#        plot = g, device = 'tiff', width =4.7,height = 2.7, compression='lzw', points= 0.8)
ggsave(sprintf('res/source_core_maas_overlap_scatter_labeled_v2.tiff'),
       # plot = g, device = 'tiff', width =6,height = 4, compression='lzw', points= 0.8)
       plot = g, device = 'tiff', width =4.7,height = 2.6, compression='lzw', points= 0.8)

# cor.test(df$coef_rural, df$coef_dx, test_type='spearman')
# 
# df$Direction = make.names(Direction)
# df2 = df[df$Direction!='Not significant\nin CD vs\nRural controls\n',]
# 
# write.table(x = df2, file = sprintf('res/rural_dx_maas_overlap_scatter_q%s_table.txt',qval_cut), quote = F, sep = '\t', row.names = F)
