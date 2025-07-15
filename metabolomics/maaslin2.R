## ** only do once - installs the required packages
# if(!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Maaslin2")
source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')
library(ggrepel)
library(Maaslin2)

# writing the necessary text file for maaslin2 run
write_maaslin2 = function(taxa, metadata_cols, maas_path)
{
  # seprating taxa to taxonomy and metadata
  pos = which(names(taxa) == 'Description')
  metadata = taxa[, c('SampleID',metadata_cols) ]
  tx = taxa[ ,c( 1,(pos+1):dim(taxa)[2] ) ]
  
  browser
  # writing to text files
  metadata_file = sprintf('%s/metadata.tsv',maas_path)
  tx_file = sprintf('%s/taxonomy.tsv',maas_path)
  write.table(metadata, metadata_file, sep = '\t',quote = F,row.names = F)
  write.table(tx, tx_file, sep = '\t',quote = F,row.names = F)
}

n_name = 'CORE'

type = 'Feces'
# type = 'Serum'

## ** metadata file path. sent you the file I use. If you want to use another one change the path here
map_file = '../../multiomics/data/CORE_screening_multiomics_map_v8.txt'
## ** biom file converted to text file fitted for R.
## ** again, sent you the files I use. If you want to use others change the path here
if (type == 'Feces')
  biom_file = '../data/feces-mbx12-s106-f417_biom.txt'
if (type == 'Serum')
  biom_file = '../data/serum-normalized-no-top-2-s114-f391_biom.txt'


map = read.table(map_file,sep="\t", header=TRUE, comment.char = '', quote = '')
biom = read.table(biom_file,sep="\t", header=TRUE, comment.char = '',skip = 1, row.names = 1, quote = '')
biom = as.data.frame(t(biom))
row.names(biom) = make.names(row.names(biom))

if (type == 'Feces')
{
  map = map[!is.na(map$ID_MBX_stool),]
  map$SampleID = make.names(map$ID_MBX_stool)
}
if (type == 'Serum')
{
  map = map[!is.na(map$ID_MBX_serum),]
  map$SampleID = make.names(map$ID_MBX_serum)
}
biom = biom[map$SampleID,]

biom$SampleID = row.names(biom)
map$SampleID = make.names(as.character(map$SampleID))
map$Description = map$SampleID

taxa = merge(map, biom, by='SampleID', all.x = T)
names(taxa) = make.names(names(taxa))

## ** directory the results will be saved to
path = 'maasin2_res/'

head_path = sprintf('%s/',path)
dir.create(head_path)

## ** name of this specific run, will create a sub directory in the main maaslin results directory 
name = sprintf('%s_%s_Dx_gender_age',n_name, type )
## ** fixed_cols are the metadata parameters toy want to check here
fixed_cols = c('Dx','Gender','Age')
## ** if a parameter you want to check is categoric and has more than 2 options, you need to decide 
## ** what is there reference you want to compare to (your base). here for example I decided 
## ** to compare Disease_Status (CD_active, CD_active, Control) to Control. 
## ** if there is more than one such parameter use ; between them.
# referce_cols = 'Disease_Status,Control'
# random_cols = c('pn_ID')
random_cols = c()

metadata_cols = c(fixed_cols, random_cols)

maas_path = sprintf('%s/%s',head_path, name)
dir.create(maas_path)

write_maaslin2(taxa, metadata_cols, maas_path)

out_path = sprintf('%s/res/',maas_path)

fit_data <- Maaslin2( input_data = sprintf('%s/taxonomy.tsv',maas_path),
                      input_metadata = sprintf('%s/metadata.tsv',maas_path),
                      output = out_path,
                      fixed_effects = fixed_cols,
                      # reference = referce_cols,
                      random_effects = random_cols)


qval_cut = 0.25

df = fit_data$results
df = df[df$metadata == 'Dx',]

col = ifelse(df$qval <= qval_cut, 'Up','NS')
col[col == 'Up' & df$coef<0] = 'Down'


lbl = df$feature
lbl[df$qval > qval_cut] = ''
volc_p = ggplot(df, aes(x=coef, y=-log10(qval), label = lbl)) + 
  geom_point(aes(colour = col)) + 
  # geom_text_repel(size = 2, seed = 4, nudge_y = 0.2, max.overlaps = 2, box.padding = 0.1) + 
  geom_text_repel(size = 2, seed = 4, nudge_y = 0.2, box.padding = 0.4) + 
  scale_color_manual(values = c('blue','gray','red'), name='') + 
  theme_bw() + theme(panel.grid = element_blank())
volc_p
ggsave(sprintf('%s/%s_q%s_volcano.tiff', out_path, 'Dx', qval_cut),
       plot = volc_p, device = 'tiff', width =5,height = 3.5, compression='lzw')

