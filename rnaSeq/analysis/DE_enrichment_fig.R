library(ggplot2)
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')


data_file = 'data/CORE_DE_toppfun_table_yael.txt'
na_str = c('no_data','_','NA','unknown', 'other','na','No_followup','')
df = read.table(data_file,sep="\t", header=TRUE, na.strings = na_str, comment.char = '', quote = '')

short_name_file = 'data/CORE_DE_toppfun_table_choose.txt'
dfc = read.table(short_name_file,sep="\t", header=TRUE, na.strings = na_str, comment.char = '', quote = '')
dfc$used_in_wgcna[dfc$short_name == 'extracellular matrix'] = F

# df$include_old = df$include
df$include = df$ID %in% dfc$ID


df = merge(df, dfc[,c('ID','short_name','used_in_wgcna')], by= 'ID', all.x = T, all.y = F)

## add the single cell we want to add
pos = df$ID == '219c5d0cde7f6082755154f54db221413ec555cb'
df$include[pos] = T
df$short_name[pos] = 'goblet cells'
df$used_in_wgcna[pos] = T


df = df[ order(df$include, decreasing = T),]
df = df[ order(df$used_in_wgcna, decreasing = T),]
df = df[df$dataset %in% c('Pathway','GO: Biological Process','GO: Cellular Component',
                          'GO: Molecular Functio','Coexpression Atlas',
                          'Mouse Phenotype','ToppCell Atlas'),]

write.table(df, 'data/CORE_DE_toppfun_table_yael_v2.txt', sep = '\t', quote = F, row.names = F)

# df = df[!df$Name %in% c('brush border membrane','inflammatory response'),]

# df$include[df$Name %in% df$Name[df$include == 'TRUE']] = 'TRUE'

df = df[df$include == 'TRUE',]
temp = as.data.frame(table(df$Name))
# df = df[df$Name %in% temp$Var[temp$Freq>1],]

df$FDR_BH = as.numeric(df$FDR_BH)

df = df[df$dataset %in% c('GO: Molecular Functio',
                          'GO: Biological Process',
                          'GO: Cellular Component',
                          'Pathway','ToppCell Atlas'),]

df$group[df$group == 'active_up'] = 'Active > Remission'
df$group[df$group == 'active_down'] = 'Remission > Active'
df$group[df$group == 'remission_up'] = 'Remission > Control'
df$group[df$group == 'remission_down'] = 'Control > Remission'

df$group = factor(df$group, levels = rev(c('Active > Remission',
                                       'Remission > Active',
                                       'Remission > Control',
                                       'Control > Remission')))
cols = rev(c('#A04747','#ff8080','#ffc080','#8d8dff'))

df$Name = sprintf('%s (%s)', df$short_name, df$ID)
df$Name[df$short_name == 'goblet cells'] = 'goblet cells'

df2 = maditr::dcast(data = df,formula = Name~group,fun.aggregate = sum,value.var = "FDR_BH")
# df2 = maditr::dcast(data = df,formula = short_name~group,fun.aggregate = sum,value.var = "FDR_BH")
df2 = as.data.frame(df2)
row.names(df2) = df2$Name
# row.names(df2) = df2$short_name
df2 = df2[,-1]
df2 = ifelse(df2==0,0,1)
ord <- hclust( dist(df2, method = "canberra"), method = "ward.D" )$order
df$Name = factor(df$Name, row.names(df2)[ord])
# df$short_name = factor(df$short_name, row.names(df2)[ord])



df$FDR_BH = -1*log10(df$FDR_BH)

y_breaks = levels(df$Name)
y_col = ifelse(y_breaks %in% df$Name[df$used_in_wgcna == T], '#8A2D3B','black')

g = ggplot(df, aes(x=group, y=Name, size =FDR_BH, colour = group)) + 
  geom_point() + theme_bw() + 
  geom_point(aes(fill = dataset), alpha = 0) + 
  # scale_color_brewer(palette = 'Set1', name = '') + 
  scale_color_manual(values = cols, name = '', guide = guide_legend(reverse = TRUE)) + 
  scale_fill_brewer(palette = 'Set3', name = '') + 
  # scale_shape_manual(values = c(16,17,18, 15), name = '') + 
  # facet_grid(dataset~., scales = 'free', space = 'free') + 
  scale_y_discrete(breaks = y_breaks) +  
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(color = y_col),
        panel.grid.major = element_blank()) + 
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 22, size = 3))) +
  # ggside::geom_ysidetile(aes(x = '', yfill = dataset, ycolour = dataset)  ) + 
  xlab('') + ylab('') + labs(size='-log10(FDR)') 
  
df2 = as.data.frame(df2)
df2$Name = row.names(df2)
df2$Name = factor(df2$Name, df2$Name[ord])
df2 = merge(df2, df[,c('Name','dataset')], by = 'Name')
leg = ggplot(df2, aes(y = Name, x = 0)) + 
  geom_point(aes(color = dataset), shape = 15, size = 3, show.legend = F) + 
  theme_classic() + 
  scale_colour_brewer(palette = 'Set3') +
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) 
  
# scale_colour_manual(values = unique_cols) #

# annotate hm with the bar
g2 = g + annotation_custom(ggplotGrob(leg), 
                          xmin = -.1, xmax = .4, 
                          ymin = 0, ymax =length(unique(df2$Name))+0.5) # n + .5 since tiles are 1 in width centered on the value
g2

# ggsave('DE_function_fig.pdf', plot = g, device = 'pdf', width = 7.3,height = 5.5)
ggsave('DE_function_fig.pdf', plot = g2, device = 'pdf', width = 6.5,height = 6)
# ggsave('DE_function_fig_full.pdf', plot = g, device = 'pdf', width = 11.5,height = 7.5)

g_flip = g+coord_flip() + 
  annotation_custom(ggplotGrob(leg + coord_flip()), 
                    xmin = -.1, xmax = .4, 
                    ymin = 0, ymax =length(unique(df2$Name))+0.5) + 
  theme(axis.text.x= element_text(hjust = 1, vjust = 0.5, color = y_col),
        axis.ticks.x = element_blank())
ggsave('DE_function_fig_flip.pdf', plot = g_flip, device = 'pdf', width =8.5,height = 4.3)

temp = vector(length = dim(df2)[1], mode = 'character')
for ( i in 1:length(temp) )
  temp[i] = sprintf('%s %s %s %s', df2[i,'Active > Remission'],
                    df2[i,'Remission > Active'],
                    df2[i,'Remission > Control'],
                    df2[i,'Control > Remission'])


# df3 = df
# df3$FDR_BH[grepl('<', df3$group)] = -1*df3$FDR_BH[grepl('<', df3$group)]
# df3$group = gsub('[<>]','vs', df3$group)
# 
# dir = ifelse(df3$FDR_BH < 0,'Up','Down')
# dir = factor(dir,c('Up','Down'))
# g_bar = ggplot(df3, aes(x=FDR_BH, y = Name, fill = dir)) + 
#   geom_bar(stat="identity", position='dodge') +
#   geom_vline(xintercept = 0, col = "gray") + 
#   # geom_col(position = position_dodge(preserve = "single"), width = 1.2) + 
#   # facet_grid(FDR_BH < 0~., scales = 'free', space = 'free') + 
#   facet_grid(.~group, scales = 'free') + 
#   scale_fill_brewer(palette = 'Set1', name = '') + 
#   theme_bw() + xlab('Direction * FDR') + ylab('') + 
#   theme(panel.grid.major.x = element_blank(), 
#         panel.grid.minor = element_blank())
# g_bar
# 
# df4 = df3
# df4$FDR_BH[df4$group == 'Active vs Remission'] = -1*df4$FDR_BH[df4$group == 'Active vs Remission'] 
# df4$group[df4$group == 'Active vs Remission'] = 'Remission vs Active'
# g_bar2 = ggplot(df4, aes(x=FDR_BH, y = Name, fill = group)) + 
#   geom_bar(stat="identity", position=position_dodge2(preserve = "single")) +
#   geom_vline(xintercept = 0, col = "gray") + 
#   # geom_col(position = position_dodge(preserve = "single"), width = 1.2) + 
#   # facet_grid(FDR_BH < 0~., scales = 'free', space = 'free') + 
#   # facet_grid(.~group, scales = 'free') + 
#   scale_fill_brewer(palette = 'Set1', name = '') + 
#   theme_bw() + xlab('Direction * FDR') + ylab('') + 
#   theme(panel.grid.major.x = element_blank(), 
#         panel.grid.minor = element_blank())
# 
# g_bar2
