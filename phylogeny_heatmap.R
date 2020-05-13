library(ggplot2)
library(ggtree)
library(data.table)
library(ape)
library(tidyverse)
library(lubridate)
library(scales)
library(treedater)

#load full tree and all dates. Perform rough dating using default treedater
ml_tree <- ape::read.tree("~/Documents/covid_deep/data/consensus/MA_combined_genomes_with_ref.fasta.treefile")
ml_tree$tip.label <- gsub("_Con.*", "", ml_tree$tip.label)
all_dates <- fread("~/Documents/covid_deep/data/consensus/combined_dates.csv", 
                   data.table = FALSE, header = FALSE, col.names = c("sample", "date"))
all_dates <- all_dates[!all_dates$date %in% c("", "None"),]
date_vector <- lubridate::decimal_date(lubridate::as_date(all_dates$date))
names(date_vector) <- all_dates$sample
samples_missing_date <- ml_tree$tip.label[!ml_tree$tip.label %in% all_dates$sample]
ml_tree <- ape::drop.tip(ml_tree, tip = samples_missing_date)
dated_tree <- dater(ml_tree, date_vector, s=29903, omega0 = 9e-3)



# tree_file <- "~/Documents/covid_deep/data/consensus/transcluster_out/index_seqs.fasta.timetree.nwk"
# tree <- ggtree::read.tree(tree_file)
# tree$tip.label <- gsub("_.*", "", gsub("cluster_", "", tree$tip.label))

cluster_file <- "~/Documents/covid_deep/data/consensus/transcluster_out/transclusters.csv"
clusters <- data.table::fread(cluster_file, data.table = FALSE) %>% as_tibble()
clusters$calendar_date <- lubridate::as_date(clusters$calendar_date)
clusters <- clusters[order(clusters$calendar_date),]
index_samples <- clusters$sample[!duplicated(clusters$cluster)]

tree <- ape::drop.tip(dated_tree, dated_tree$tip.label[!dated_tree$tip.label %in% index_samples])
tree$tip.label <- clusters$cluster[match(tree$tip.label, clusters$sample)]
class(tree) <- "phylo"

# collect cluster counts
clusters$calendar_date <- lubridate::ceiling_date(clusters$calendar_date, "day")
cluster_counts <- clusters %>% group_by(cluster, calendar_date) %>%
  summarise(count = n())

# can round to weeks ext with 
ncluster <- length(unique(clusters$cluster))
all_dates <- lubridate::ceiling_date(seq(min(clusters$calendar_date), max(clusters$calendar_date), 1), "day")
all_dates <- all_dates[!duplicated(all_dates)]
count_matrix <- matrix(0, nrow = ncluster, ncol = length(all_dates), 
                       dimnames = list(unique(clusters$cluster), all_dates))
count_matrix[cbind(cluster_counts$cluster, as.character(as.numeric(lubridate::ceiling_date(cluster_counts$calendar_date, "day"))))] <- cluster_counts$count

# generate phylogenetic heatmap
p <- ggtree(tree, mrsd=min(all_dates)) + 
  geom_tiplab(size=2, align=TRUE, linesize=.5) + 
  theme_tree2()

cluster_counts$dd <- as.numeric(cluster_counts$calendar_date)

date_breaks <- unique(lubridate::ceiling_date(cluster_counts$calendar_date, "week"))
date_breaks <- date_breaks[order(date_breaks)]

p + geom_facet(panel = "Transmission Cluster Counts", data = cluster_counts, geom = geom_tile, 
               mapping=aes(x = dd, fill = count), size=2) + 
  scale_fill_fermenter(palette="RdYlBu", breaks=c(1,2,5,10,50,100)) + 
  theme(axis.text.x= element_text(size=12, angle=45, hjust = 1),
        strip.text.x = element_text(size = 12)) +
  scale_x_continuous(breaks = as.numeric(date_breaks), 
                     labels =as.character(date_breaks))

