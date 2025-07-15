source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(ggplot2)
library(ggrepel)

# name = 'active_vs_remission'
# DE_file = 'res/CORE_S_plus_geneFiltered_CD_Active_(n=27)_vs_CD_remission_(n=55)_deseq2_res.tsv'

name = 'remission_vs_control'
DE_file = 'res/CORE_S_plus_geneFiltered_CD_remission_(n=55)_vs_Control_(n=36)_deseq2_res.tsv'

df = read.table(DE_file, header = T, sep = '\t', comment.char = ,quote = '')

FDR_cutoff = 0.05
FC_cutoff = log2(1.5)
col = ifelse( abs(df$log2FoldChange) >= FC_cutoff & df$padj <= FDR_cutoff, 'Up','NS')
col[col == 'Up' & df$log2FoldChange<0] = 'Down'

# cols = c('#e31a4c','gray','#4a65c5')
# if (name == 'remission_vs_control')
#   cols = c('#fbba3a','gray','#4a65c5')
# if (name == 'active_vs_remission')
#   cols = c('#e31a4c','gray','#fbba3a')

cols = c('#ff8080','gray','#8d8dff')
if (name == 'remission_vs_control')
  cols = c('#ffc080','gray','#8d8dff')
if (name == 'active_vs_remission')
  cols = c('#ff8080','gray','#ffc080')

lbl = gsub('_.*','',df$Gene)
lbl[df$padj > FDR_cutoff] = NA
col = factor(col,levels = c('Up','NS','Down'))
volc_p = ggplot(df, aes(x=log2FoldChange, y=-log10(padj), label = lbl)) + 
# volc_p = ggplot(df, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(aes(colour = col), size=1) + 
  # geom_text_repel(size = 2, seed = 4, nudge_y = 0.2, max.overlaps = 2, box.padding = 0.1) + 
  # geom_text_repel(size = 2, seed = 4, max.overlaps = 2, box.padding = 0.1) + 
  # scale_color_manual(values = c('blue','gray','red'), name='') + 
  scale_color_manual(values = cols, name='') + 
  ylab('-log10(FDR)') + xlab('log2(FC)') + 
  theme_bw() + theme(panel.grid = element_blank())
volc_p
out_path = 'res/'
ggsave(sprintf('%s/%s_q%s_FC%s_volcano.tiff', out_path, name, FDR_cutoff, 'log2(1.5)'),
       plot = volc_p, device = 'tiff', width =3.5,height = 2, compression='lzw')

