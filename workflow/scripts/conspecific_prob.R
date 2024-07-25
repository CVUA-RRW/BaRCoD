require(tidyverse)

# Inputs
distance_table <- snakemake@input[["distance_table"]]
min_id <- snakemake@params[["id"]]

# Outputs
output <- snakemake@output[["conspecific_prob_table"]]
graph <- snakemake@output[["conspecific_plot"]]

min_id <- as.numeric(min_id)*100
print(id)


tab <- 
    read_tsv(distance_table) %>%
    mutate(same_taxid = if_else(query_taxid == target_taxid, 1,0)) %>%
    group_by(id_cut = cut(id, breaks=seq(min_id,100, 0.5))) %>%
    summarize(n_pairs = n(),
              n_same_tax = sum(same_taxid),
              conspecific_prob = n_same_tax/n_pairs) 

tab %>%
    write_tsv(output)

tab %>% 
    ggplot(aes(x=id_cut, y=conspecific_prob))+
    geom_point()+
    geom_line()+
    xlab("Identity interval") +
    ylab("Conspecific probability") 

ggsave(graph)