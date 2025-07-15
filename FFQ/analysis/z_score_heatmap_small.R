source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

metadata_prms = c('Cohort','Dx','Dx_activity','new_diagnosis',
                  'Age','Gender','BMI','Smoking_YN',
                  'pn_ID')
col_parm = 'Dx_activity'

add_name = ''
# add_name = sprintf('%s_maaslin_sig', add_name)

ffq_file = '../data/CORE_FFQ_with_metadata.txt'
ffq = read.table(ffq_file, header = T, sep = '\t', row.names = 1)

# ## check samples
# map_file = '../../multiomics/data/CORE_screening_multiomics_map_v4.txt'
# metadata_df = read.table(map_file, header = T, sep = '\t', na.strings = c('NA','','#N/A'), comment.char = '',quote = '')
# metadata_df = metadata_df[!is.na(metadata_df$ID_FFQ),]
# all(sort(row.names(ffq)) == sort(metadata_df$pn_ID))

food_prms = names(ffq)[!names(ffq) %in% metadata_prms]

maas_file = '../res/maaslin2/Dx_activity_age_gender_ctlRef/res/significant_results.tsv'
maas_res = read.table(maas_file, header = T, sep = '\t',)
maas_res = maas_res[maas_res$metadata == 'Dx_activity',]
if (grepl('maaslin_sig',add_name))
  food_prms = unique(maas_res$feature)

cor_type = 'spearman'
groups = unique(ffq[[col_parm]])

ffq = ffq[ffq[[col_parm]] %in% groups, ]

out_path = '../res/z_heatmaps/'
dir.create(out_path)

df2 = ffq[food_prms]
# row.names(df2) = df$pn_ID

good_pos = c()
for ( i in 1:dim(df2)[2] )
{
  if ( sum(!is.na(df2[,i])) > 5 & sum(!df2[,i] ==0, na.rm = T )!=0 )
  {
    good_pos = c(good_pos, i)
  }
}
df2 = df2[,good_pos]

good_pos = c()
for ( i in 1:dim(df2)[1] )
{
  if ( sum(!is.na(df2[i,])) > 2  )
  {
    good_pos = c(good_pos, i)
  }
}
df2 = df2[good_pos,]


df3 = df2
for ( i in 1:dim(df3)[2] ) ## z score
{
  # df3[,i] = log10(df3[,i]) + 0.0001
  df3[,i] = (df3[,i] - mean(df3[,i], na.rm = T)) / sd(df3[,i], na.rm = T)
  # df3[,i] = (df3[,i] / mean(df3[,i]))
}
# df3 = log(df3+0.0001)

cols20 = c('#f3730c', '#114efe', '#ce1c20', '#4363d8', '#f58231', '#911eb4',
           '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
           '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
           '#000075', '#808080', '#ffffff', '#000000')
dx_active_cols= c('#fbba3a', '#4a65c5', '#e31a4c')

Group = ffq[[col_parm]]
gs = unique(Group)
for( i in 1:length(gs))
  Group[Group == gs[i]] = dx_active_cols[i]



# Gender = ffq$Dx
# Gender[Gender == 'Control'] = '#344CB7'
# Gender[Gender == 'CD'] = '#8E1616'
# Gender[Gender == 'UC'] = '#FFD65A'

Gender = ffq$Gender
Gender[Gender == 'Male'] = '#81BFDA'
Gender[Gender == 'Female'] = '#FF90BC'

col2 = colorRampPalette(rev( c('#67001F', '#B2182B', '#D6604D', '#F4A582', '#FDDBC7',
                               '#FFFFFF', '#D1E5F0', '#92C5DE', '#4393C3', '#2166AC', '#053061') ))(200)

#Load latest version of heatmap.3 function
library('gplots')
library(devtools)
source_url('https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R')

mydist=function(c) {dist(c,method='euclidian')}
myclust=function(c) {hclust(c,method='ward.D')}
temp = t(as.matrix(df3))

temp[temp > abs(min(temp)) ] = abs(min(temp))

# width = 7
# height =  15
width = 6.3
height = 4.7
if (grepl('maaslin_sig',add_name))
  height=5.7
# temp[temp > 10] = 10
# pdf(file=sprintf('%s/z_score_heatmap_%s.pdf', out_path, paste(groups, collapse = '_') ), width = 7, height = 15)
# width=20
pdf(file=sprintf('%s/z_score_kcal_norm_heatmap_%s%s.pdf', 
                 out_path, paste(groups, collapse = '_'), add_name ), 
    width = width, height = height)
# pdf(file=sprintf('%s/z_score_absolute_heatmap_%s.pdf', out_path, paste(groups, collapse = '_') ), width = 7, height = 15)
par(cex.main=1)
heatmap.3(temp,  col = col2, ColSideColors  = cbind(Group, Gender),
          distfun = mydist, hclustfun = myclust,
          na.rm = TRUE, scale='none', dendrogram='both', 
          margins=c(6,14),
          # margins=c(3,3),
          Rowv=TRUE, Colv=TRUE, symbreaks=F, density.info='none', trace='none',
          # ain= sprintf('z_score_heatmap_%s', paste(groups, collapse = '_')), labCol=F,
          cexRow=1, 
          ColSideColorsSize=1.4, side.height.fraction = 0.5 )

##  heatmap(temp,col = col2)
# legend(x=0, y=0.8,legend=c(gs,
#                            '',
#                            'Female','Male'),
#        fill=c(cols20[1:length(gs)],
#               NA,
#               '#F37199','#66D2CE'), 
#        # border=F, 
#        bty='n', y.intersp = 0.7, cex=0.7)
dev.off()

# pdf(file='heatmap3_example.pdf')
# main_title='Drug Response Predictions'
# par(cex.main=1)
# heatmap.3(prob_matrix, hclustfun=myclust, distfun=mydist, na.rm = TRUE,
#           scale='none', dendrogram='both', margins=c(6,12),
#           Rowv=TRUE, Colv=TRUE, ColSideColors=clab, RowSideColors=rlab,
#           symbreaks=FALSE, key=TRUE, symkey=FALSE,
#           density.info='none', trace='none', main=main_title, labCol=FALSE,
#           labRow=drug_names, cexRow=1, col=rev(heat.colors(75)),
#           ColSideColorsSize=7, RowSideColorsSize=2, KeyValueName='Prob. Response')
# legend('topright',legend=c('Basal','LumA','LumB','Her2','Claudin','Normal','','Positive','Negative','NA','','Targeted','Chemo','','Approved','Experimental'),
#        fill=c('red','blue','cyan','pink','yellow','green','white','black','white',
# 'grey','white','darkorchid','darkred','white','green','darkgreen'), border=FALSE, bty='n', y.intersp = 0.7, cex=0.7)
# dev.off()


# heatmap3(t(as.matrix(df3)),  col = col2, ColSideColors  = cbind(Group, colGender), balanceColor = F,
#          method = 'ward.D',# distfun = 'euclidean',
#          legend('topright',legend=c('Basal','LumA','LumB'),
#                 fill=c('red','blue','cyan')))

# heatmap(t(as.matrix(df3)),  col = col2, ColSideColors  = Group)
# 
# legend('topright',legend=c('Basal','LumA','LumB'),
#                 fill=c('red','blue','cyan'))



