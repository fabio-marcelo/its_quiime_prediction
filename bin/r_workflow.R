# install
q = c("BiocManager")
for (pacote in q) {
  if(!require(q,character.only = TRUE)) utils::install.packages(q)
  library(q,character.only = TRUE)
}



requiredPackages = c("UpSetR", "fastqcr", "pander", "dplyr",
                     "tidyr", "kableExtra", "remotes", "quarto", "data.table",
                     "labeling", "ggplot2", "ggpubr", "rmarkdown", "R.utils",
                     "gridExtra", "plyr", "cowplot", "vegan",
                     "tidyverse", "scales", "maditr","ape",
                     "devtools","ggdendro","gridExtra","knitr",
                     "pander","plotly","png","tidyverse","vegan",
                     "tidyverse", "picante", "Polychrome",
                     "RColorBrewer", "scales")   


for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p,dependencies = T)
  library(p,character.only = TRUE)
}


Packages = c("ngsReports", "phyloseq", "DESeq2", "ShortRead", "Biostrings",
             "ComplexHeatmap", "microbiome", "dada2", "DECIPHER")                      
for(o in Packages){
  #if(!require(o,character.only = TRUE, quietly = TRUE)) BiocManager::install(o)
  library(o,character.only = TRUE)
}

devtools::install_github("jbisanz/qiime2R")
devtools::install_github('https://github.com/smped/ngsReports/tree/RELEASE_3_17')
remotes::install_github('https://github.com/YuLab-SMU/ggtree',dependencies = T)
remotes::install_github('https://github.com/yiluheihei/microbiomeMarker/releases/tag/v0.0.1')
BiocManager::install('karyoploteR')

library(c(qiime2R, ngsReports, ggtree, microbiomeMarker, karyoploteR))

# clear workspace
rm(list = ls())

caminho <- '/home/mapa/nextflow_run/fastqc_dir'
files_raw <- list.files(caminho,
                        pattern = "fastqc.zip$", full.names = T)


fdl_raw <- ngsReports::FastqcDataList(files_raw)
summary_raw <- ngsReports::getModule(fdl_raw, "Summary")
reads_raw <- ngsReports::readTotals(fdl_raw)

save.image(file='/home/mapa/nextflow_run/r_workflow.RData')