library(ggplot2)
source("/pita/users/tzipi/code/R_figs/figs_funcs.R")

core_screening_file = '../data/CORE_screening_multiomics_map_v8.txt'
map = read.table(file = core_screening_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'), quote = '', comment.char = '')


out_path = 'plots/Dx_boxplots/'

indexes = c('faith_pd','dysbiosis_index',
            'NOVA_class','I.MEDAS_SCORE',
            'Single.cell.atlas..goblet.cells','GO..UDP.glycosyltransferase',
            'REACTO..mucins.glycosylation','GO..immune.adaptive..T.cells.',
            'REACTO..immune.innate..granulocyte.','GO..cytokine.activity',
            'GO..anibacterial..epithelia.'# ,'GO..extracellular_matrix'
            )

group1 = NA
group2 = NA

test_name = 'Dx_activity'


cols = c( '#e31a4c', '#fbba3a', '#4a65c5')
for ( index in indexes )
{
  map_f = map[!is.na(map[[index]]), ]
  for (d in unique(map_f[[test_name]]))
  {
    map_f[[test_name]][map_f[[test_name]] == d] = 
      sprintf('%s\n(n=%s)',d, sum(map_f[[test_name]] == d, na.rm = T))
  }
  g = ggplot(map_f) + 
    # geom_violin(aes_string(test_name, index), fill = 'blue', color = NA, alpha = 0.8) + 
    geom_boxplot(aes_string(test_name, index, fill  = test_name)) + 
    scale_fill_manual(values = cols, name = '') + 
    theme_bw() + xlab('') + 
    theme(panel.grid = element_blank(), legend.position = 'none')
  res = add_significant_asterix_to_plot_BH_v2(p = g, 
                                              var = as.factor(map_f[[test_name]]), val = map_f[[index]], 
                                              test_type = 'wilcox')
  # res2 = ggaddsig(p = g, var = map_f[[test_name]],
  #                 val = map_f[[index]])
  g = res[[1]]
  
  out_file = sprintf('%s/%s_%s_boxplot.tiff', out_path, test_name, index)
  ggsave(out_file, g ,device = 'tiff', width = 4,height = 4, compression = 'lzw')
}

nova_mocus_scatter = ggplot(map, aes(NOVA_class, 
       REACTO..mucins.glycosylation, colour = Dx_activity)) + 
  geom_point() + 
  ylab('Mucins glycosylation') + 
  xlab('NOVA class') + 
  geom_smooth(method = 'lm', se = F) + 
  scale_colour_manual(values = cols, name = '') + 
  theme_bw() + 
  theme(panel.grid = element_blank())

ggsave('plots/nova_mocus_scatter.tiff', nova_mocus_scatter,
       device = 'tiff', width = 4.5,height = 2.5, compression = 'lzw')

medas_mocus_scatter = ggplot(map, aes(I.MEDAS_SCORE, 
                                     REACTO..mucins.glycosylation, 
                                     colour = Dx_activity)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F) + 
  ylab('Mucins glycosylation') + 
  xlab('I-MEDAS score') + 
  scale_colour_manual(values = cols, name = '') + 
  theme_bw() + 
  theme(panel.grid = element_blank())

ggsave('plots/medas_mocus_scatter.tiff', medas_mocus_scatter,
       device = 'tiff', width = 4.5,height = 2.5, compression = 'lzw')

