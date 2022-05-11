require(tidyverse)
require(scales)

# Inputs
distances <- snakemake@input[["distance_table"]]
consensus <- snakemake@input[["cons"]]

# Outputs
dist_dist <- snakemake@output[["dist_plot"]]
cons_dist <- snakemake@output[["consensus_plot"]]

theme_set(theme_bw())

dist_tbl <- read_tsv(distances) %>%
    ggplot(aes(x=distance)) +
        geom_bar(position = "dodge2") +
        xlab("Pairwise distance")

ggsave(dist_dist)

cons_tbl <- read_tsv(consensus) %>%
    pivot_longer(ends_with('mismatches'), 
                 names_to="mismatches",
                 names_pattern="([0-9]+) mismatches",
                 values_to="rank") %>%
    ggplot(aes(x=mismatches, fill=rank)) +
        geom_bar() 
        
ggsave(cons_dist)