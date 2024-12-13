#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)

# args[1] <- "~/workspace/2022/gapsnake.test/genomes/"

args[1] <- gsub("/$","",args[1])

if(args[1] == "")
  stop("No valid directory with input genomes provided.")

if(!dir.exists(args[1]))
  stop(paste("Directory",args[1],"does not exist."))

valid_extensions <- "\\.fasta$|\\.fasta\\.gz$|\\.fna$|\\.fna\\.gz$|\\.faa$|\\.faa\\.gz$|\\.fa$|\\.fa\\.gz$|\\.fn$|\\.fn\\.gz$"

genomes <- dir(args[1], pattern = valid_extensions, recursive = TRUE)

dt <- data.table(sample = basename(gsub(valid_extensions,"", genomes)),
                 genome_file = paste0(args[1],"/",genomes),
                 taxonomy = "auto",
                 biomass = "auto",
                 translation_table = 11)
dt[, medium := paste0("models/",sample,"/",sample,"-medium.csv")]

if(any(duplicated(dt$sample)))
  stop("Sample names not unique. There are genome files in different sub-directories with the same basename...")

fwrite(dt, "samples.tsv", sep = "\t", quote = FALSE)
