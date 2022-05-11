require(tidyverse)
require(scales)


# Inputs
primers <-  snakemake@input[["primers"]]
barcodes <- snakemake@input[["barcodes"]]

# Outputs
primer_hist <- snakemake@output[["primer_hist_plot"]]
barcode_hist <- snakemake@output[["barcode_hist_plot"]]
size_hist <- snakemake@output[["size_hist_plot"]]

# Read data
primer_tbl <- read_tsv(primers)
barcode_tbl <- read_tsv(barcodes, col_names=c('seqid', 'taxid', 'start', 'end', 'length'))


theme_set(theme_bw())

# Primer histogram - number of binding site per primer
primer_tbl %>%
    mutate(primer=str_split_fixed(query, "\\[", 2)[,1]) %>%
    ggplot(aes(x=primer)) +
        geom_bar(position = "dodge2")

ggsave(primer_hist)


# Barcode Histogram - number of barcode per taxid
barcode_tbl %>%
    group_by(taxid) %>%
    summarize(count=n()) %>%
    ggplot(aes(x=count)) +
        geom_bar(position = "dodge2") +
        xlab("Barcode per taxid")

ggsave(barcode_hist)


# Size Histogram - sequence length distribution
barcode_tbl %>% 
    ggplot(aes(x=length))+
        geom_bar(position="dodge2") +
        xlab("Barcode length")

ggsave(size_hist)