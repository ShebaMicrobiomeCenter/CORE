source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')

taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-29_merge/res/16S_DB1-29_merged/taxonomy/taxonomy_gg2_source_risk.tsv'

asv_file = '../data/core_asvs_list.txt'
df = read.table(asv_file, header = T)

res = sv_2_tax_v2(svs = df$ASV, taxonomy_file = taxonomy_file)

write.table(x = res, file = '../data/core_asvs_list_taxa.txt', quote = F, sep = '\t', row.names = F)