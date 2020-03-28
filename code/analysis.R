require(data.table)
require(preprocessCore)
require(qsmooth)
require(tidyverse)
require(BuenColors)

read_counts <- function(filename){
  file <- read_tsv(filename)[1:2]
  colnames(file) <- c('ensg', 'expected_count')
  file 
}
files <- list.files(".", pattern="ReadsPerGene.out.tab")
l <- lapply(files, read_counts)
fnames <- files %>% gsub("ReadsPerGene.out.tab","",. ) %>%
                    gsub("_1","-1",.) %>%
                    gsub("_2","-2",.) %>% 
                    gsub("_3","-3",.)
names(l) <- fnames

df <- l %>% 
  bind_rows(., .id="id") %>%
  spread(id, expected_count)

reads_per <- df %>% dplyr::select(-ensg) %>%
  colSums(.) %>%
  data.frame(n_reads = .) %>%
  rownames_to_column("cell")

p1 <- ggplot(reads_per, aes(x = cell, y = n_reads)) +
  geom_bar(stat = "identity")

qs <- qsmooth(as.matrix(df[,-1]), fnames)

norm_df <- qsmoothData(qs) %>%
  as.tibble(.) %>%
  mutate(ensg = df$ensg)

reads_per2 <- norm_df %>% dplyr::select(-ensg) %>%
  colSums(.) %>%
  data.frame(n_reads = .) %>%
  rownames_to_column("cell")

p2 <- ggplot(reads_per2, aes(x = cell, y = n_reads)) +
  geom_bar(stat = "identity")

combined_cells <- norm_df %>%
  gather(id, expected_count, -ensg) %>%
  mutate(cell_type = strsplit(id, split="-") %>%
           lapply(., function(x) x[[1]]) %>%
           unlist(.)) %>%
  group_by(cell_type, ensg) %>%
  summarise(count = sum(expected_count)) %>%
  spread(cell_type, count)

combined_cells = subset(combined_cells, select = -c(hs_BFU,hs_CD34,hs_CFU))
names(combined_cells) <- c('ensg','CB_BFUE','CB_CD34','CB_CFUE',"CB_eBaso","CB_lBaso","CB_ortho","PB_BFUE", 'PB_CD34','PB_CFUE','PB_eBaso','PB_lBaso','PB_ortho','PB_poly','PB_pro','CB_poly','CB_pro')

combined_cells <- separate(data = combined_cells, col = ensg, into = c("ensg", "remove"), sep = "\\.")
combined_cells <- subset(combined_cells,select=-c(remove))
# cpm 
combined_cpm <- combined_cells %>%
  mutate_if(is.numeric, cpm, log = F)
#log2
log2_cpm_p1 <- function(x) log2(x + 1)
combined_cpm_l2 <- combined_cpm %>%
  mutate_if(is.numeric, log2_cpm_p1)
# make both cpm long
combined_cpm_long <- combined_cpm %>%
  gather(cell_types, cpm, -ensg)
combined_cpm_l2_long <- combined_cpm_l2 %>%
  gather(cell_types, cpm, -ensg)
# plot function
plot_gene <- function(long_df, gene){
  gene <- long_df %>% filter(ensg == gene)
  print(gene)
  p1 <- ggplot(gene, aes(x = cell_types, y = cpm)) +
    geom_bar(stat = "identity")
  print(p1)
}

# Test

plot_gene(combined_cpm_long, "ENSG00000107890")
plot_gene(combined_cpm_long, "ENSG00000095787")
plot_gene(combined_cpm_l2_long, "ENSG00000107890")
plot_gene(combined_cpm_l2_long, "ENSG00000095787")
# final plot looks like log2 cpm will be best
THE_PALETTE <- jdb_palette("solar_rojos")
gene <- combined_cpm_l2_long %>% filter(ensg %in% c("ENSG00000107890", "ENSG00000095787")) %>%
  mutate(ensg = factor(ensg))
ggplot(gene, aes(x = cell_types, y = cpm, fill = cpm)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  facet_wrap(~factor(ensg), ncol = 2) +
  pretty_plot() + 
  scale_fill_gradientn(colors = THE_PALETTE)

saveRDS(combined_cells, "count.rds")

PostScriptTrace('trial3.eps', 'trial3.xml',
                charpath=TRUE, charpos=FALSE,
                setflat=NULL, defaultcol="black",
                encoding="ISO-8859-1", scaleEPS=.01)