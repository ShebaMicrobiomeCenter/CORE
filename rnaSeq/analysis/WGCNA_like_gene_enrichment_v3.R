library(WGCNA)
library(ggplot2)
options(stringsAsFactors = FALSE);

source('/pita/users/tzipi/code/R_figs/WGCNA_funcs.R')

modules_order = c('GO: immune adaptive (T cells)',
                  'REACTOME: immune innate (granulocyte)',
                  'GO: cytokine activity',
                  'GO: anibacterial (epithelia)',
                  'GO: extracellular_matrix','Single cell atlas: goblet cells',
                  'GO: UDP glycosyltransferase','REACTOME: mucins glycosylation')

# add_omic = 'Metabolites'
add_omic = 'FFQ'
add_name = add_omic

# index_type = 'eigengene'
index_type = 'PC1'
# index_type = 'zscore_mean'

# min_q = 0.1
min_q = 0.25

# cor_type = 'spearman'
cor_type = 'WGCNA_default'


out_path = 'res/WGCNA_like_enrichemnt/'
dir.create(out_path)
out_path = sprintf('%s/CORE_S_plus/', out_path)
dir.create(out_path)

# in_path = '/pita/users/tzipi/projects/multiomics/SOURCE/multi_analysis/WGCNA/rnaSeq_main/res_v2/'
# israel_source_name = sprintf('SOURCE_Israel_TI_withFFQv10_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_%s','WGCNA_default')
# israel_source_short_name = 'SOURCE_Israel_TI'
# israel_ti = load( sprintf('%s/%s/data.RData',in_path, israel_source_name) )
# module_df = df


## load CORE side data
tpm_file = 'res/CORE_S_plus_geneFiltered_txi_res_geneFiltered.txt'
ftr_df = read.table(tpm_file, header = T, sep = '\t', na.strings = c('NA',''), row.names = 1)
ftr_df = as.data.frame(t(ftr_df))
names(ftr_df) = gsub('__','_',names(ftr_df))

## setting enrichment genes instead of modules gene list
enrichment_genes_file = 'data/core_DE_enrichment_genes.txt'
eg = read.table(enrichment_genes_file, header = T, sep = '\t', comment.char = '',quote = '')
names(eg)[1] == 'Gene'
temp = data.frame(Gene = names(ftr_df), short_gene = gsub('_.*','',names(ftr_df)))
temp = temp[temp$short_gene %in% eg$Original_Symbol,]
eg = eg[eg$Original_Symbol %in% temp$short_gene,]
eg = merge(eg, temp, by.x = 'Original_Symbol', by.y = 'short_gene', all.x = T, all.y = F)
# eg$module_color = sprintf('%s: %s', eg$direction, eg$term_short)
eg$module_color = sprintf('%s', eg$term_short)
module_df = eg 

# remove unwanted terms 
module_df = module_df[!module_df$term %in% c('abnormal_self_tolerance;_MP:0005005',
                                             'abnormal_granulocyte_physiology;_MP:0002462',
                                             'extracellular_matrix;_GO:0031012'),]

name = sprintf('CORE_S_plus_by_%s', 'enrichment_genes')
map_file = '../multiomics/data/CORE_screening_multiomics_map_v8.txt'
metadata_df_basic = read.table(map_file, header = T, sep = '\t', na.strings = c('NA','','#N/A'), comment.char = '',quote = '')
metadata_df_basic = metadata_df_basic[!is.na(metadata_df_basic$ID_rnaSeq_TI),]
row.names(metadata_df_basic) = make.names(metadata_df_basic$ID_rnaSeq_TI)

## load org WGCNA metadata, from source("set_WGCNA_metadata.R")
metadata_df = read.table('data/WGCNA_map.txt', header = T, sep = '\t', na.strings = c('NA','','#N/A'), comment.char = '',quote = '')
metadata_cat_df = read.table('data/WGCNA_map_cats.txt', header = T, sep = '\t', na.strings = c('NA','','#N/A'), comment.char = '',quote = '')

# metadata_df = metadata_df[ metadata_df$Active_vs_control != 1 | is.na(metadata_df$Active_vs_control),]
# add_name = sprintf('%s_noActive',add_name)

# wanted_cats = c('General','Bacteria', add_omic)
wanted_cats = c('General', add_omic)
temp = make.names( metadata_cat_df$var[metadata_cat_df$cat %in% wanted_cats ] )
metadata_df = metadata_df[, temp ]
add_vars = metadata_cat_df$var[metadata_cat_df$cat == add_omic ]

row.names(ftr_df) = make.names(row.names(ftr_df))
row.names(metadata_df) = make.names(row.names(metadata_df))
ftr_df = ftr_df[row.names(metadata_df),]

metadata_prms = names(metadata_df)

## sync genes between the datasets
good_genes = intersect(names(ftr_df), module_df$Gene)
module_df = module_df[module_df$Gene %in% good_genes, ]
ftr_df = ftr_df[, module_df$Gene]

## ME calculation
# Calculate MEs with color labels
# WGCNA_plot_genes_correlaitons_by_module(ftr_df, module_df$module_color, out_path)

if (index_type == 'eigengene')
{
  res = WGCNA_calculate_module_eigengenes_and_corrs(ftr_df,module_df$module_color, metadata_df, cor_type, name)
} else if (index_type == 'zscore_mean')
{
  res = WGCNA_calculate_module_zscore_mean_and_corrs(ftr_df,module_df$module_color, metadata_df, cor_type, name)
  name = sprintf('%s_%s',name, index_type)
} else if (index_type == 'PC1')
{
  res = WGCNA_calculate_module_PC1_and_corrs(ftr_df,module_df$module_color, metadata_df, cor_type, name)
  name = sprintf('%s_%s',name, index_type)
}
MEs = res[[1]]
moduleTraitCor = res[[2]]
moduleTraitPvalue = res[[3]]

me = MEs
names(me) = gsub('ME','',names(me))
write.table(me, sprintf('%s/%s_module_MEs.txt',out_path, name), sep = '\t', row.names = T, quote = F)

# # my version
# heatmap = WGCNA_ME_met_heatmap(moduleTraitCor, moduleTraitPvalue, name, metadata_prms)
# # ggsave(sprintf('%s/%s_module_cor_heatmap%s.pdf',out_path, name, add_name),plot = heatmap, device = 'pdf', width = 3.5+dim(moduleTraitPvalue)[1]/2.5,height = 3+dim(moduleTraitPvalue)[2]/10)



t = metadata_prms[!make.names(metadata_prms) %in% make.names(add_vars)]
heatmap_s_fdr = WGCNA_ME_met_heatmap_q(moduleTraitCor,
                                       moduleTraitPvalue, name,
                                       metadata_prms = t,
                                       wanted_modules_order = modules_order,
                                       cluster_cols_order_flag = T,
                                       min_q = min_q)
# heatmap_s_fdr = WGCNA_ME_met_heatmap_pq(moduleTraitCor, 
#                                        moduleTraitPvalue, name, 
#                                        metadata_prms = t,
#                                        wanted_modules_order = modules_order, 
#                                        cluster_cols_order_flag = T, 
#                                        min_q = min_q, min_p = 0.05)
# res = WGCNA_ME_met_heatmap_full_res(moduleTraitCor, moduleTraitPvalue, name, metadata_prms)

heatmap_s = heatmap_s_fdr[[1]]
y_lables = ggplot_build(heatmap_s)$layout$panel_params[[1]]$y$get_labels()

ggsave(sprintf('%s/%s_module_cor_heatmap_%s_q%s.pdf',out_path, name, add_name, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10, limitsize = F)
# # ggsave(sprintf('%s/%s_module_cor_heatmap_p%s_q%s_spearman.pdf',out_path, name, min_p, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10, limitsize = F)

heatmap_s_table = heatmap_s_fdr[[2]]
write.table(heatmap_s_table, sprintf('%s/%s_module_cor_heatmap_%s_q%s.txt',out_path, name, add_name, min_q), sep = '\t', row.names = F, quote = F)


mdc = heatmap_s_table
col2 = colorRampPalette(rev( c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                               "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") ))(200)

temp = levels(mdc$variable)
ffq_vars = temp[temp %in% add_vars]
mdc$variable = as.character(mdc$variable)
mdc$variable[mdc$variable == 'Smoking_YN'] = 'Smoking'
mdc$variable[mdc$variable == 'city_socioeconomi_rank_2019'] = 'Socioeconomi rank'
mdc$variable[mdc$variable == 'FCP'] = 'Fecal calprotectin^'
mdc$variable[mdc$variable == 'CRP'] = 'CRP^'
mdc$variable[mdc$variable == 'dysbiosis_index_stool'] = 'Dysbiosis index'
mdc$variable[mdc$variable == 'faith_pd_stool'] = 'Faith\'s PD'
mdc$variable[mdc$variable == 'Active_vs_remission'] = 'CD active vs remission'
mdc$variable[mdc$variable == 'Remission_vs_control'] = 'CD remission vs control'
mdc = mdc[mdc$variable!='Dx',]

clinical_vars = c('CD remission vs control','CD active vs remission',
                  'CRP^','Fecal calprotectin^')
dem_vars = c('Age','Gender','BMI','Smoking','Socioeconomi rank')
taxa_vars = c('Dysbiosis index','Faith\'s PD')

y_vars = c(clinical_vars, dem_vars, taxa_vars, ffq_vars)
mdc$variable = factor(mdc$variable, levels = rev(y_vars))

mdc$category = NA
mdc$category[mdc$variable %in% clinical_vars] = 'Clinical'
mdc$category[mdc$variable %in% dem_vars] = 'Demo-\ngraphics'
mdc$category[mdc$variable %in% taxa_vars] = 'Microbial'
mdc$category[mdc$variable %in% ffq_vars] = 'FFQ'

mdc$category = factor(mdc$category, levels = c('Clinical','Demo-\ngraphics',
                                               'Microbial','FFQ'))
heatmap = ggplot(mdc, aes(x=module,y=variable,fill = correlation)) +
  geom_tile() + 
  geom_tile(data = mdc[mdc$q_value <= min_q & 
                         mdc$variable %in% ffq_vars,], 
            color = 'black', size = 0.2, aes(x=module,y=variable,fill = correlation)) +
  scale_fill_gradientn(colours = col2, limits = c(-1,1)) +
  # geom_point(aes(shape = p_value<=0.05)) + scale_shape_manual(values = c(NA, 8)) +
  # geom_text(aes(label=sprintf('%.2f\n%.2e',correlation,p_value )), size=3) + 
  scale_colour_manual(values = c('black','white')) + 
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) + 
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) + coord_flip() +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  ylab('') + 
  geom_text(aes(label=signif(p_value, 1)), size=2.5, color = 'gray20') + 
  theme_classic() +
  facet_grid(category~., scales = 'free', space = 'free') + 
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.y = element_text(angle = 0),
        panel.spacing = unit(0.1,'lines') ,
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.3) ,
        strip.background = element_rect(colour = 'black', fill = 'gray90')
        ) 

heatmap

ggsave(sprintf('%s/%s_module_cor_heatmap_%s_q%s_v2.pdf',out_path, name, add_name, min_q),
       plot = heatmap, device = 'pdf', width = 4+dim(moduleTraitPvalue)[1]/3,height = 2.5+length(y_lables)/10, limitsize = F)

# ggsave(sprintf('%s/%s_heatmap_modules%s_%s_FDR_small.pdf',out_path, name, add_name, cor_type),
#        plot = g, device = 'pdf', width = 3,height = 3,limitsize = FALSE )


## check specific corrs
# sanity check
all(row.names(metadata_df)== row.names(MEs))
pos = !is.na(metadata_df$FCP) & !is.na(metadata_df$NOVA_class)

cor.test(metadata_df$NOVA_class[pos], metadata_df$FCP[pos], method = 'spearman')
cor.test(metadata_df$NOVA_class[pos], MEs$`MEREACTOME: mucins glycosylation`[pos], method = 'spearman')
cor.test(MEs$`MEREACTOME: mucins glycosylation`[pos], metadata_df$FCP[pos], method = 'spearman')

cor.test(metadata_df$NOVA_class[pos], metadata_df$FCP[pos])
cor.test(metadata_df$NOVA_class[pos], MEs$`MEREACTOME: mucins glycosylation`[pos])
cor.test(MEs$`MEREACTOME: mucins glycosylation`[pos], metadata_df$FCP[pos])


cor.test(metadata_df$I.MEDAS_SCORE[pos], metadata_df$FCP[pos])
cor.test(metadata_df$I.MEDAS_SCORE[pos], MEs$`MEREACTOME: mucins glycosylation`[pos])
cor.test(MEs$`MEREACTOME: mucins glycosylation`[pos], metadata_df$FCP[pos])


cor.test(metadata_df$NOVA_class[pos], metadata_df_old$dysbiosis_index_stool[pos])
cor.test(metadata_df$NOVA_class[pos], MEs$`MEREACTOME: mucins glycosylation`[pos])
cor.test(MEs$`MEREACTOME: mucins glycosylation`[pos], metadata_df_old$dysbiosis_index_stool[pos])


cor.test(metadata_df$I.MEDAS_SCORE[pos], metadata_df_old$dysbiosis_index_stool[pos])

