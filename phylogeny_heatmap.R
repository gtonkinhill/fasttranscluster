library(ggplot2)
library(ggtree)
library(data.table)
library(ape)
library(tidyverse)
library(lubridate)
library(scales)

tree <- ape::read.tree("./test_out/index_seqs.fasta.treefile")
clusters <- data.table::fread("./test_out/transclusters.csv", data.table = FALSE) %>% as_tibble()
# clusters$cluster <- clusters$cluster + 1
clusters$calendar_date <- lubridate::as_date(clusters$calendar_date)
interval_size <- 1

# root tree using earliest sequence (it might be better to use an outgroup)
# node_dates <- clusters$calendar_date[match(gsub(".*_", "", tree$tip.label), as.character(clusters$sample))]
# root_node <- tree$tip.label[which.min(node_dates)]
# tree <- ape::root(tree, root_node)

# relabel tree tips with cluster
tree$tip.label <- gsub("_.*", "", gsub("cluster_", "", tree$tip.label))

# collect cluster counts
clusters$calendar_date <- lubridate::ceiling_date(clusters$calendar_date, "week")
cluster_counts <- clusters %>% group_by(cluster, calendar_date) %>%
  summarise(count = n())

# can round to weeks ext with 
ncluster <- length(unique(clusters$cluster))
all_dates <- lubridate::ceiling_date(seq(min(clusters$calendar_date), max(clusters$calendar_date), 1), "week")
all_dates <- all_dates[!duplicated(all_dates)]
count_matrix <- matrix(0, nrow = ncluster, ncol = length(all_dates), 
                       dimnames = list(unique(clusters$cluster), all_dates))
count_matrix[cbind(cluster_counts$cluster, as.character(as.numeric(lubridate::ceiling_date(cluster_counts$calendar_date, "week"))))] <- cluster_counts$count

# generate phylogenetic heatmap
p <- ggtree(tree) + 
  geom_tiplab(size=2, align=TRUE, linesize=.5) + 
  theme_tree2()

cluster_counts$dd <- as.numeric(cluster_counts$calendar_date)

date_breaks <- unique(lubridate::ceiling_date(cluster_counts$calendar_date, "month"))
date_breaks <- date_breaks[order(date_breaks)]

p + geom_facet(panel = "Transmission Cluster Counts", data = cluster_counts, geom = geom_point, 
               mapping=aes(x = dd, color = count), size=2) + 
  scale_colour_fermenter(palette="RdYlBu") + 
  theme(axis.text.x= element_text(size=12, angle=45, hjust = 1),
        strip.text.x = element_text(size = 12)) +
  scale_x_continuous(breaks = as.numeric(date_breaks), 
                     labels =as.character(date_breaks))
