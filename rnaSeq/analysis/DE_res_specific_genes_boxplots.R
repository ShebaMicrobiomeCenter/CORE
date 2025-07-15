source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(ggplot2)
library(patchwork)

ctl = 'Control'
rem = 'CD_remission'
act = 'CD_Active'

# 'active_vs_remission'
DE_file_ar = 'res/CORE_S_plus_geneFiltered_CD_Active_(n=27)_vs_CD_remission_(n=55)_deseq2_res.tsv'
de_ar = read.table(DE_file_ar, header = T, sep = '\t', comment.char = ,quote = '')
de_ar$Gene_short = gsub('_.*','',de_ar$Gene)

# 'active_vs_control'
DE_file_ac = 'res/CORE_S_plus_geneFiltered_CD_Active_(n=27)_vs_Control_(n=36)_deseq2_res.tsv'
de_ac = read.table(DE_file_ac, header = T, sep = '\t', comment.char = ,quote = '')
de_ac$Gene_short = gsub('_.*','',de_ac$Gene)

# 'remission_vs_control'
DE_file_rc = 'res/CORE_S_plus_geneFiltered_CD_remission_(n=55)_vs_Control_(n=36)_deseq2_res.tsv'
de_rc = read.table(DE_file_rc, header = T, sep = '\t', comment.char = ,quote = '')
de_rc$Gene_short = gsub('_.*','',de_rc$Gene)

metadata_file = '../multiomics/data/CORE_screening_multiomics_map_v8.txt' 
map = read.table(metadata_file, header = T, sep = '\t', comment.char = '',quote = '')
map = map[!is.na(map$ID_rnaSeq_TI),]
row.names(map) = make.names(map$ID_rnaSeq_TI)

outcome_file = '../multiomics/data/CORE_outcome.txt'
odf = read.table(file = outcome_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'), quote = '', comment.char = '')
map = merge(map, odf, by='pn_ID', all.x = T, all.y = F)
row.names(map) = make.names(map$ID_rnaSeq_TI)

tpm_file = 'res/CORE_S_plus_geneFiltered_txi_res_geneFiltered.txt'
tpm = read.table(tpm_file, header = T, sep = '\t', na.strings = c('NA',''), row.names = 1)
tpm = as.data.frame(t(tpm))

map = map[row.names(tpm),]

group = 'Immune adaptive\n(T cells)'
# group = 'Immune innate\n(granulocytes/myeloid)'
# group = 'Anti-bacterial\n(epithelia)'
# group = 'Mucus glycosylation\n(epithelia)'
# group = 'Cytokine activity'


if (group == 'Immune adaptive\n(T cells)')
{
  group_md = 'GO..immune.adaptive..T.cells.'
  group_md_name = 'GO: Immune adaptive\n(T cells)'
  # wanted_genes = c('TLR4','STAT4','NLRP1','LILRB4')
  wanted_genes = c('STAT4','ITK','TNFSF8','DOCK8')
}
if (group == 'Immune innate\n(granulocytes/myeloid)')
{
  group_md = 'REACTO..immune.innate..granulocyte.'
  group_md_name = 'REACTO: Immune\ninnate granulocyte'
  wanted_genes = c('S100A8','S100A9','OSM','TREM1')
}
if (group == 'Anti-bacterial\n(epithelia)')
{
  group_md = 'GO..anibacterial..epithelia.'
  group_md_name = 'GO: Anti-bacterial\nepithelia'
  wanted_genes = c('DUOX2','CEACAM6','DMBT1','REG1B')
}
if (group == 'Mucus glycosylation\n(epithelia)')
{
  group_md = 'REACTO..mucins.glycosylation'
  group_md_name = 'REACTO: Mucins\nglycosylation'
  wanted_genes = c('MUC1','B3GNT7','GALNT4','GALNT7')
}
if (group == 'Cytokine activity')
{
  group_md = 'GO..cytokine.activity'
  group_md_name = 'GO: Cytokine\nactivity'
  wanted_genes = c('AATF','AATF','AATF','AATF')
}

set_wanted_gene_plot = function(wanted_gene, tpm, map, de_ar, de_rc, act, ren, ctl)
{
  tpm_short_gene = gsub('__.*','',names(tpm))
  
  pos = which(tpm_short_gene == wanted_gene)
  ## handle 0
  val = log2( tpm[[pos]] + (min(tpm[[pos]][!tpm[[pos]]==0], na.rm = T)/2)  )
  g = ggplot(tpm) + 
    scale_fill_manual(values = c("#ff8080", "#ffc080", "#8d8dff"), name = '') + 
    geom_boxplot(aes(y=map$Dx_activity, x=val, fill = map$Dx_activity), 
                 outlier.size = 0.7) +
    theme_bw() + xlab('log2(TPM)') + ylab(wanted_gene) + 
    theme(panel.grid = element_blank(), legend.position = 'none',
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
  ac_q = de_ac$padj[de_ac$Gene_short == wanted_gene]
  ar_q = de_ar$padj[de_ar$Gene_short == wanted_gene]
  rc_q = de_rc$padj[de_rc$Gene_short == wanted_gene]

  manual_p = data.frame(var1 = c(act, rem, act), 
                        var2 = c(rem, ctl, ctl),
                        p=c(ar_q, rc_q, ac_q))
  res = ggaddsig(p=g, var = map$Dx_activity, 
                 val = val, 
                 manual_pval = manual_p, 
                 angle = 90, label_size = 5, scale_var = 5, line_width = 0.5,
                 asterix_scale_var = 1.5, next_step_scale = 3, )
  return(res[[1]])
  
}


p1 = set_wanted_gene_plot(wanted_genes[1], tpm, map, de_ar, de_rc, act, ren, ctl)
p2 = set_wanted_gene_plot(wanted_genes[2], tpm, map, de_ar, de_rc, act, ren, ctl)
p3 = set_wanted_gene_plot(wanted_genes[3], tpm, map, de_ar, de_rc, act, ren, ctl)
p4 = set_wanted_gene_plot(wanted_genes[4], tpm, map, de_ar, de_rc, act, ren, ctl)

val = 
p5 = ggplot(map) + 
  scale_fill_manual(values =  c("#ff8080", "#ffc080", "#8d8dff"), name = '') + 
  geom_boxplot(aes_string(y='Dx_activity', x=group_md, fill = 'Dx_activity'), 
               outlier.size = 0.7) +
  theme_bw() + xlab(group_md_name) + ylab('') + 
  theme(panel.grid = element_blank(), legend.position = 'none',
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
p5_res = ggaddsig(p=p5, var = map$Dx_activity, 
                        val = map[[group_md]], 
                        angle = 90, label_size = 5, scale_var = 5, line_width = 0.5,
                        asterix_scale_var = 1.5, next_step_scale = 3, )
p5 = p5_res[[1]]

p = (p1 + ggtitle(group) + theme(axis.title.x = element_blank())) + 
  (p2 + theme(axis.title.x = element_blank())) + 
  (p3 + theme(axis.title.x = element_blank())) + 
  p4 + 
  (ggplot() + theme_minimal()) + 
  p5 + 
  plot_layout(ncol = 1)

ggsave(sprintf('plots/terms_genes_boxplots/%s_boxplots.pdf', 
               gsub('/','_',group, fixed = T)), 
       p, device = 'pdf', width =2.1,height = 7.2)


