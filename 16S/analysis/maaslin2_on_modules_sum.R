source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(ggplot2)

out_path = '../res_stool_TI_v3_amnonClean//R_res/maaslin2_modules/'

# maas_path_string = '../res_stool_TI_v3_amnonClean//R_res/maaslin2_modules//CORE_stool_age_gender_Dx_activity_%s/res/all_results.tsv'
# out_file = sprintf('%s/CORE_stool_age_gender_Dx_activity_sum_bar.pdf',out_path )
# out_file2 = sprintf('%s/CORE_stool_age_gender_Dx_activity_heatmap.pdf',out_path )

maas_path_string = '../res_stool_TI_v3_amnonClean//R_res/maaslin2_modules//CORE_stool_age_gender_%s/res/all_results.tsv'
out_file = sprintf('%s/CORE_stool_age_gender_sum_bar.pdf',out_path )
out_file2 = sprintf('%s/CORE_stool_age_gender_heatmap.pdf',out_path )
out_table = sprintf('%s/CORE_stool_age_gender_heatmap_table.txt',out_path )

# maas_path_string = '../res_stool_TI_v3_amnonClean//R_res/maaslin2_modules//CORE_stool_CD_age_gender_%s/res/all_results.tsv'
# out_file = sprintf('%s/CORE_stool_CD_age_gender_sum_bar.pdf',out_path )
# out_file2 = sprintf('%s/CORE_stool_CD_age_gender_heatmap.pdf',out_path )

# maas_path_string = '../res_stool_TI_v3_amnonClean//R_res/maaslin2_modules//CORE_stool_noActive_age_gender_%s/res/all_results.tsv'
# out_file = sprintf('%s/CORE_stool_noActive_age_gender_sum_bar.pdf',out_path )
# out_file2 = sprintf('%s/CORE_stool_noActive_age_gender_heatmap.pdf',out_path )

# maas_path_string = '../res_stool_TI_v3_amnonClean//R_res/maaslin2_modules//CORE_stool_rem_age_gender_%s/res/all_results.tsv'
# out_file = sprintf('%s/CORE_stool_rem_age_gender_sum_bar.pdf',out_path )
# out_file2 = sprintf('%s/CORE_stool_rem_age_gender_heatmap.pdf',out_path )

modules_order = c('GO: immune adaptive (T cells)',
                  'REACTO: immune innate (granulocyte)',
                  'GO: cytokine activity',
                  'GO: anibacterial (epithelia)',
                  # 'GO: extracellular_matrix',
                  'Single cell atlas: goblet cells',
                  'GO: UDP glycosyltransferase','REACTO: mucins glycosylation')

module_names = make.names(modules_order)

sum_df = data.frame(module = module_names, sig_n = NA, sig_up_n = NA, sig_dn_n = NA)

maas_df = data.frame()
for (md in module_names)
{
  maas_path = sprintf(maas_path_string, md)
  df = read.table(maas_path, header = T, sep = '\t')
  df = df[df$metadata == md, ]
  maas_df = rbind(maas_df, df)
  df = df[df$qval <= 0.25, ]
  sum_df$sig_n[sum_df$module == md] = dim(df)[1]
  sum_df$sig_up_n[sum_df$module == md] = sum(df$coef > 0)
  sum_df$sig_dn_n[sum_df$module == md] = sum(df$coef < 0)
}

for (i in 1:length(modules_order))
{
  sum_df$module[sum_df$module == make.names(modules_order[i])] = modules_order[i] 
  maas_df$metadata[maas_df$metadata == make.names(modules_order[i])] = modules_order[i] 
}

g = ggplot(sum_df, aes(y=module, x=sig_n)) + 
  geom_bar(stat = "identity") + 
  theme_bw() # + theme(axis.text.x = element_text(angle = 45, hjust = 1))

sum_df2 = reshape::melt(sum_df, id.vars = 'module', variable_name = 'direction')
sum_df2$module = factor(sum_df2$module, modules_order)
# g2 = ggplot(sum_df2, aes(y=module, x=value, fill = direction)) + 
#   geom_bar(position="dodge", stat="identity") + 
#   theme_bw()
# g2


g3 = ggplot(sum_df2[sum_df2$direction!='sig_n',], 
            aes(x=module, y=value, fill = direction)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_bw() + ylab('# significant maaslin results') + xlab('') + 
  scale_fill_manual(values = c('#E55050','#3d8dbf'), labels = c('Positive','Negative'),name = 'Correlation') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank()) 
g3


ggsave(out_file, plot = g3, device = 'pdf', width = 2.8,height = 4.6)


col2 = colorRampPalette(rev( c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                              "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") ))(200)

wanted_features = unique(maas_df$feature[maas_df$qval<=0.25])
maas_df = maas_df[maas_df$feature %in% wanted_features, ]
maas_df$metadata = factor(maas_df$metadata, rev(modules_order))

maas_df$feature = gsub('.sp[0-9]+_','_',maas_df$feature)
maas_df$feature = gsub('_[0-9]+','',maas_df$feature)
temp = maditr::dcast(data = maas_df,formula = metadata~feature,fun.aggregate = sum,value.var = "coef")
temp = as.data.frame(temp)
row.names(temp) = temp$metadata
temp = temp[,names(temp) != 'metadata']
ord <- hclust( dist(t(temp), method = "euclidean"), method = "ward.D" )$order
maas_df$feature = factor(maas_df$feature, levels = names(temp)[ord])

max_coef = max(abs(maas_df$coef))

gm = ggplot(maas_df, aes(y=metadata,x=feature,fill = coef)) +
 geom_tile() +
 scale_fill_gradientn(colours = col2, limits = c(-1*max_coef,max_coef)) +
 # geom_point(aes(shape = p_value<=0.05)) + scale_shape_manual(values = c(NA, 8)) +
 # geom_text(aes(label=sprintf('%.2f\n%.2e',correlation,p_value )), size=3) +
 scale_colour_manual(values = c('black','white')) +
 # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) +
 theme(axis.text.x = element_text(angle = 45, hjust=1, size =8),
       legend.box = "horizontal") +
 # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) + coord_flip() +
 scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
 geom_point(aes(y=metadata,x=feature, alpha = qval <= 0.25), shape = 8, size = 1) +
 scale_alpha_manual(values = c(0,1), name= 'Q<=0.25') +
 xlab('') + ylab('') # + guides(alpha = "none")
 # theme(axis.text.y=element_blank(), 
 #       axis.ticks.y=element_blank()) 
# ggsave(out_file2, plot = gm, device = 'pdf', width = 7.5,height = 3+length(wanted_features)/10)
ggsave(out_file2, plot = gm, device = 'pdf', width = 6+length(wanted_features)/8,height = 4)


write.table(x = maas_df, file = out_table, sep = '\t', row.names = F, quote = F)
