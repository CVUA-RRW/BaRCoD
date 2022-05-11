require(tidyverse)
require(scales)

# Inputs
clusters <- snakemake@input[["sizes"]]

# Outputs
size_dist <- snakemake@output[["size_plot"]]
taxid_size <- snakemake@output[["taxid_size"]]

theme_set(theme_bw())

size_tbl <- read_tsv(clusters) %>%
    ggplot(aes(x=cluster_size)) +
        geom_bar(position = "dodge2")

ggsave(size_dist)


size_tbl <- read_tsv(clusters) %>%
    group_by(taxid) %>%
    summarize(count=n()) %>%
        ggplot(aes(x=count)) +
        geom_bar(position = "dodge2") +
        xlab("Unique barcodes per taxid")

ggsave(taxid_size)