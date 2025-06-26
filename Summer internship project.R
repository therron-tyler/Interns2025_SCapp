library(shiny)
library(ggplot2)
expr_data <- read.delim("Sarcoidosis_MacSubset_metadata_BatchInfo.tsv",
                        header = TRUE, row.names = 1, check.names = FALSE)
meta_data <- read.delim("Sarcoidosis_MacSubset_metadata_BatchInfo.tsv",
                        header = TRUE, row.names = 1, check.names = FALSE)
head(expr_data)
head(meta_data)
sum(is.na(expr_data))
sum(is.na(meta_data))


