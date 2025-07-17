library(ggplot2)
library(patchwork)

core_screening_file = '../data/CORE_screening_multiomics_map_v8.txt'
map = read.table(file = core_screening_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'), quote = '', comment.char = '')

terms_index_type = 'eigengene'
# terms_index_type = 'zscore_mean'
if (terms_index_type == 'zscore_mean')
{
  me_file = '../../rnaSeq/res/WGCNA_like_enrichemnt/CORE_S_plus/CORE_S_plus_by_enrichment_genes_zscore_mean_module_MEs.txt'
  me = read.table(file = me_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'), quote = '', comment.char = '')
  add_vars = names(me)
  map = map[,!names(map) %in% add_vars]
  me$SampleID = make.names(row.names(me))
  
  me_f = me[,c('SampleID',add_vars)]
  map$ID_rnaSeq_TI[!is.na(map$ID_rnaSeq_TI)] = make.names(map$ID_rnaSeq_TI[!is.na(map$ID_rnaSeq_TI)])
  map = merge(map,me_f, by.x = 'ID_rnaSeq_TI', by.y = 'SampleID', all.x = T, all.y = F )
}

out_path = 'plots/subgroup_boxplots_v3/'

indexes = c(# 'faith_pd','dysbiosis_index',
            # 'NOVA_class','I.MEDAS_SCORE',
            'GO..immune.adaptive..T.cells.','REACTO..immune.innate..granulocyte.',
            'GO..cytokine.activity','GO..anibacterial..epithelia.',
            'Single.cell.atlas..goblet.cells',
            'GO..UDP.glycosyltransferase','REACTO..mucins.glycosylation')
indexes_name = c(# 'faith_pd','dysbiosis_index',
  # 'NOVA_class','I.MEDAS_SCORE',
  'GO: Immune adaptive\n(T cells)','REACTO: Immune innate\n(granulocyte)',
  'GO: Cytokine\nactivity','GO: Anibacterial\n(epithelia)',
  'Single cell atlas:\nGoblet cells',
  'GO: UDP\nglycosyltransferase','REACTO: Mucins\nglycosylation')



group1 = NA
group2 = NA

# test_name = 'CD_remission_vs_CD_Active'
test_name = 'biologics_vs_treatment_naive'
# test_name = 'remission_biologics_vs_no_or_5asa'
# test_name = 'remission_LS_under_350_vs_LS_over_350'
# test_name = 'remission_endoscopic_non_inflammad_vs_endoscopic_inflammad'

if (test_name == 'CD_remission_vs_CD_Active')
{
  group1 = map$Dx_activity == 'CD_remission' & !is.na(map$Dx_activity)
  group2 = map$Dx_activity == 'CD_Active'
  cols = c("#ff8080", "#ffc080", "#8d8dff")
}
if (test_name == 'biologics_vs_treatment_naive')
{
  group1 = map$any_bioligics_CD == 1 & !is.na(map$any_bioligics_CD)
  group2 = map$new_diagnosis == 'Diagnosis'
  cols = c("#ff8080", "#ffc080", "#8d8dff")
}
if (test_name == 'remission_biologics_vs_no_or_5asa')
{
  group1 = map$any_bioligics_remission == 1 & !is.na(map$any_bioligics_remission) 
  group2 = map$any_bioligics_remission == 0 & !is.na(map$any_bioligics_remission) 
  cols = c("#FF9B45", "#FED16A", "#8d8dff")
}
if (test_name == 'remission_endoscopic_non_inflammad_vs_endoscopic_inflammad')
{
  group1 = map$endoscopic_inflammation == 0 & !is.na(map$endoscopic_inflammation) 
  group2 = map$endoscopic_inflammation == 1 & !is.na(map$endoscopic_inflammation) 
  cols = c("#FFAAAA", "#FFDCDC", "#8d8dff")
}
if (test_name == 'remission_LS_under_350_vs_LS_over_350')
{
  remission_pos = map$Dx_activity == 'CD_remission'
  group1 = map$gastro_highest_lewis_score < 350 & 
    !is.na(map$gastro_highest_lewis_score) &
    remission_pos
  
  group2 = map$gastro_highest_lewis_score > 350 & 
    !is.na(map$gastro_highest_lewis_score) &
    remission_pos
  cols = c("#FED16A", "#FFF4A4", "#8d8dff")
}
group3 = map$Dx == 'Control'
g3 = 'Control'

map[[test_name]] = NA
# map[[test_name]][group1] = sprintf('%s\n(n=%s)',
#                                    gsub('_vs_.*','',gsub('remission_','',test_name)),
#                                    sum(group1) )
# map[[test_name]][group2] = sprintf('%s\n(n=%s)',
#                                    gsub('.*_vs_','',test_name),
#                                    sum(group2) )

g1 = gsub('_vs_.*','',gsub('^remission_','',test_name))
g2 = gsub('.*_vs_','',test_name)
g1 = stringr::str_to_sentence(gsub('_',' ',g1))
g2 = stringr::str_to_sentence(gsub('_',' ',g2))
g1 = gsub('^Ls','LS', g1)
g2 = gsub('^Ls','LS', g2)

map[[test_name]][group1] = g1
map[[test_name]][group2] = g2
map[[test_name]][group3] = g3

set_index_subgroup_fig = function(map, index, test_name, g1, g2, g3, cols, index_name)
{
  map_f = map[!is.na(map[[test_name]]) & !is.na(map[[index]]), ]
  
  g1_name = sprintf('%s\n(n=%s)\n',g1, sum(map_f[[test_name]] == g1) )
  map_f[[test_name]][map_f[[test_name]] == g1] = g1_name
  g2_name = sprintf('%s\n(n=%s)\n',g2, sum(map_f[[test_name]] == g2) )
  map_f[[test_name]][map_f[[test_name]] == g2] = g2_name
  g3_name = sprintf('%s\n(n=%s)\n',g3, sum(map_f[[test_name]] == g3) )
  map_f[[test_name]][map_f[[test_name]] == g3] = g3_name
  
  map_f[[test_name]] = factor(map_f[[test_name]], levels = rev(c(g3_name,g1_name,g2_name)))
  
  g = ggplot(map_f) + 
    scale_fill_manual(values = cols, name = '', guide = guide_legend(reverse = TRUE)) + 
    geom_boxplot(aes_string(y=test_name, x=index, fill = test_name), outlier.size = 0.7) +
    theme_bw() + xlab(index_name) + ylab('') + 
    theme(panel.grid = element_blank(), legend.position = 'none',
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
  res = ggaddsig(p=g, var = map_f[[test_name]], 
                 val = map_f[[index]], 
                 angle = 90, label_size = 5, scale_var = 5, line_width = 0.5,
                 asterix_scale_var = 1.5, next_step_scale = 3, )
  g=res[[1]]
  return(g)
}

i=1
p1 = set_index_subgroup_fig(map, indexes[i], test_name, g1, g2, g3, cols, indexes_name[i])
i=2
p2 = set_index_subgroup_fig(map, indexes[i], test_name, g1, g2, g3, cols, indexes_name[i])
i=3
p3 = set_index_subgroup_fig(map, indexes[i], test_name, g1, g2, g3, cols, indexes_name[i])
i=4
p4 = set_index_subgroup_fig(map, indexes[i], test_name, g1, g2, g3, cols, indexes_name[i])
i=5
p5 = set_index_subgroup_fig(map, indexes[i], test_name, g1, g2, g3, cols, indexes_name[i])
i=6
p6 = set_index_subgroup_fig(map, indexes[i], test_name, g1, g2, g3, cols, indexes_name[i])
i=7
p7 = set_index_subgroup_fig(map, indexes[i], test_name, g1, g2, g3, cols, indexes_name[i])

legend_grob <- cowplot::get_legend(p7 + theme(legend.position = 'right'))
p = (p1) + 
  (p2) + 
  (p3) + 
  (p4)+ 
  (p5) + 
  (p6) +
  (p7)+ 
  legend_grob + 
  plot_layout(nrow = 2)
  
out_file = sprintf('%s/%s_terms_boxplot_%s.tiff', out_path, test_name, terms_index_type)
ggsave(out_file, p ,device = 'tiff', width = 8,height = 4, compression = 'lzw')


pn = set_index_subgroup_fig(map, 'NOVA_class', test_name, 
                            g1, g2, g3, cols, 'NOVA class')
pn = pn + theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = 'none')+ 
  ggtitle('Diet')

pm1 = set_index_subgroup_fig(map, 'dysbiosis_index', test_name, 
                            g1, g2, g3, cols, 'Dysbiosis index')
pm1 = pm1 + theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = 'none')+ 
  ggtitle('Microbiome')

pm2 = set_index_subgroup_fig(map, 'faith_pd', test_name, 
                             g1, g2, g3, cols, 'Faith\'s PD')
pdm = pn + pm1 + pm2 + plot_layout(nrow = 1)

out_file2 = sprintf('%s/%s_diet_taxa_boxplot.tiff', out_path, test_name)
ggsave(out_file2, pdm ,device = 'tiff', width = 8,height = 2, compression = 'lzw')


