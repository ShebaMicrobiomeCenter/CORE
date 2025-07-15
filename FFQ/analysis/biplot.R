source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(ggside)

# col_parms = c('Cohort','Dx','Dx_cohort','first_sample','ffq_origin','pn_ID','pn2')
col_parms = c('Cohort','Dx','Dx_activity','new_diagnosis',
                'Age','Gender','BMI','Smoking_YN','Age_q',
              'pn_ID')
# col_parm = 'Dx_activity'
col_parm = 'Gender'
# col_parm = 'Age_q'

side_type = 'density'
# side_type = 'boxplot'

ffq_file = '../data/CORE_FFQ_with_metadata.txt'
ffq = read.table(ffq_file, header = T, sep = '\t', row.names = 1)

ffq$Age_q = 'Q3'
ffq$Age_q[ ffq$Age < quantile(ffq$Age)[2] ] = 'Q1'
ffq$Age_q[ ffq$Age > quantile(ffq$Age)[2]  & ffq$Age < quantile(ffq$Age)[3]] = 'Q2'
ffq$Age_q[ ffq$Age > quantile(ffq$Age)[4] ] = 'Q4'


ffq$pn2 = ffq$pn_ID
temp = as.data.frame(table(ffq$pn_ID))

food_prms = names(ffq)[!names(ffq) %in% col_parms]

groups = unique(ffq[[col_parm]])
cor_type = 'spearman'
# cor_type = 'pearson'

ffq = ffq[ffq[[col_parm]] %in% c(groups) & !is.na(ffq[[col_parm]]),]


out_path = 'res/'
dir.create(out_path)

df2 = ffq[,food_prms]
# row.names(df2) = df$pn_ID

good_pos = c()
for ( i in 1:dim(df2)[2] )
{
  if ( sum(!is.na(df2[,i]) ) > 5 & sum(!df2[,i] ==0, na.rm = T )!=0 )
    # if ( sum(is.na(df2[,i]) ) == 0  & sum(!df2[,i] ==0, na.rm = T )!=0 )
  {
    good_pos = c(good_pos, i)
  }
}
df2 = df2[,good_pos]

good_pos = c()
for ( i in 1:dim(df2)[1] )
{
  if ( sum(!is.na(df2[i,]), na.rm = T) != 0  )
    # if ( sum(is.na(df2[i,])) == 0  )
  {
    good_pos = c(good_pos, i)
  }
}
df2 = df2[good_pos,]

df2 = df2[ rowSums( is.na(df2) ) == 0, ]

ffq = ffq[row.names(df2),]

df3 = df2
for ( i in 1:dim(df3)[2] )
{
  df3[,i] = (df3[,i] - mean( df3[,i], na.rm=T ) ) / sd(df3[,i], na.rm = T)
  # df3[,i] = (df3[,i] / mean(df3[,i]))
}

unique_colors = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
df.pca = prcomp(df3, center = T, scale = T)
imp = summary(df.pca)$importance[2,]
data = as.data.frame(df.pca$x)
# col_parm = 'Patient_group'

# for(col_parm in col_parms)
# {

# df$Patient_group2 = factor(df$Patient_group2, levels = c('Rural_health_<50%_in_city','Rural_health_>50%_in_city','Urban_health','Chinese_crohns'))
g = ggplot(data, aes(x=PC1, y=PC2)) +
  # geom_point(colour='gray') +
  # geom_point(aes(colour = df$Patient_group, shape = df$X16d._How_long_have_you_stayed_in_a_city_in_the_last_year_.)) +
  geom_point(aes(colour = ffq[[col_parm]] )) + guides(colour=guide_legend(title=col_parm)) +
  # scale_colour_manual(values =c('#3cb44b','#ffe119','#f58231','#9b2226')) +
  xlab(sprintf('PC1 (%.2f%%)', imp[1]*100)) + ylab(sprintf('PC2 (%.2f%%)', imp[2]*100)) +
  # geom_point() +
  # geom_text(aes(label=row.names(data)), size=3)+
  # geom_text(aes(label=names(tpm2), colour = map$loc), size=3) +
  theme_bw() # + scale_colour_brewer(palette = 'Paired')
g
# ggsave(sprintf('%s/%s_%s_pca.tiff', out_path, paste(groups, collapse = '_'), col_parm ),plot = g, device = 'tiff', width = 5,height =2.5, compression  = 'lzw')

# library(devtools)
# install_github("vqv/ggbiplot")
library(ggbiplot)
gb <- ggbiplot(df.pca, obs.scale = 1, var.scale = 1, 
               groups = ffq[[col_parm]], ellipse = F, 
               circle = F, choices = c(1,2)) + #, labels = df$pn_ID)
  # scale_colour_manual(values =c('#3cb44b','#ffe119','#f58231','#9b2226'))
  theme_bw() + guides(colour=guide_legend(title=col_parm)) 

PC = df.pca
col_data = ffq[[col_parm]]
x="PC1"
y="PC2"
# PCbiplot <- function(PC, x="PC1", y="PC2", col_data = NA) 

# PC being a prcomp object
data <- data.frame(obsnames=row.names(PC$x), PC$x)
datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
mult <- min(
  (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
  (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x]))) )
datapc <- transform(datapc,
                    v1 = mult * (get(x)),
                    v2 = mult * (get(y)) )
datapc$sums = abs(datapc$v1) + abs(datapc$v2)
datapc = datapc[order(datapc$sums, decreasing = T),]
datapc = datapc[1:10,]

cols = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', 'gray90', '#808080', '#ffffff', '#000000')
if ( col_parm  == 'Dx_activity')
  cols = c( '#e31a4c', '#fbba3a', '#4a65c5')
if ( col_parm  == 'Gender')
  cols = c( '#FF90BC','#74B2E2' )
plot <- ggplot(data, aes_string(x=x, y=y)) + 
  # geom_text(alpha=.4, size=3, aes(label=obsnames))
  # geom_point(size=3, aes(colour = col_data, alpha = col_data == 'single_sample')) + theme_bw() + 
  geom_point(size=2, aes(colour = col_data), alpha = 1 ) + 
  theme_bw() + theme(panel.grid = element_blank(), 
           axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) + 
  # scale_alpha_manual(values = c(0.8,0.3)) + 
  scale_color_manual(values =cols) + 
  scale_fill_manual(values =cols) + 
  # scale_color_brewer(palette = 'Set1') + 
  guides(colour='none', fill=guide_legend(title=NULL)) + 
  xlab(sprintf('PC1 (%.2f%%)', imp[1]*100)) + ylab(sprintf('PC2 (%.2f%%)', imp[2]*100)) #  + 
  # geom_xsidedensity(aes(fill = col_data, col = col_data, alpha = 0.5))#  +
  # geom_ysidedensity(aes(fill = col_data, col = col_data, alpha = 0.5)) +
  # theme(ggside.panel.scale = 0.25)
# plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)

# library('geomtextpath')
# plot = plot + geom_textsegment(data=datapc,aes(xend = v1, yend = v2, label = varnames)) 
# plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color=muted('blue'))
library(ggrepel)
set.seed(42)
plot = plot + geom_text_repel(data=datapc, aes(x=v1, y=v2, label=varnames), size = 3, color='black', box.padding = .3)
plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color='black')
# ggsave(sprintf('%s/%s_%s_biplot_top10.tiff', out_path, paste(groups, collapse = '_'), col_parm ),
# #        plot = plot, device = 'tiff', width =6,height = 4, compression  = 'lzw')
#       plot = plot, device = 'tiff', width =8,height = 4, compression  = 'lzw')


if (side_type == 'density')
{
  plot = plot+ geom_xsidedensity(aes(fill = col_data, col = col_data), alpha = 0.5) +
    geom_ysidedensity(aes(fill = col_data, col = col_data), alpha = 0.5) +
    theme(ggside.panel.scale = 0.25)
} else if (side_type == 'boxplot')
{
  plot = plot + geom_xsideboxplot(orientation = "y", aes(col = col_data)) +
    geom_ysideboxplot(orientation = "x", aes(col = col_data)) +
    theme(ggside.panel.scale = 0.25)
}
ggsave(sprintf('../%s/biplot/FFQ_%s_biplot.tiff', out_path, col_parm ),plot = plot, device = 'tiff', width =7,height = 4, compression  = 'lzw')
# ggsave(sprintf('%s/%s_%s_%s_biplot.tiff', out_path, paste(groups, collapse = '_'), col_parm, type ),plot = gb, device = 'tiff', width =10,height = 6, compression  = 'lzw')
# fviz_pca_biplot(df.pca)
# }


g_pc1 =  ggplot(data, ) + 
  geom_boxplot(aes(y=PC1, x=col_data, fill = col_data)) + 
  theme_bw() + 
  scale_fill_manual(values =cols) + 
  theme(panel.grid = element_blank(), 
        legend.position = 'none',
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) + 
  ylab(sprintf('PC1 (%.2f%%)', imp[1]*100)) + xlab('')
res_pc1 = add_significant_asterix_to_plot_BH_v2(p = g_pc1, val = data$PC1, var = as.factor(col_data), test_type = 'wilcox')
g_pc1 = res_pc1[[1]]
ggsave(sprintf('../%s/biplot/FFQ_%s_PC1.tiff', out_path, col_parm ),
       plot = g_pc1, device = 'tiff', width =1.5,height = 3, compression  = 'lzw')

g_pc2 = ggplot(data, ) + 
  geom_boxplot(aes(y=PC2, x=col_data, fill = col_data)) + 
  theme_bw() + 
  scale_fill_manual(values =cols) + 
  theme(panel.grid = element_blank(), 
        legend.position = 'none',
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) + 
  ylab(sprintf('PC2 (%.2f%%)', imp[2]*100)) + xlab('')
res_pc2 = add_significant_asterix_to_plot_BH_v2(p = g_pc2, val = data$PC2, var = as.factor(col_data), test_type = 'wilcox')
g_pc2 = res_pc2[[1]]
ggsave(sprintf('../%s/biplot/FFQ_%s_PC2.tiff', out_path, col_parm ),
       plot = g_pc2, device = 'tiff', width =1.5,height = 3, compression  = 'lzw')

