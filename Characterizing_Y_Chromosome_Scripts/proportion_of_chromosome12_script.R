library(tidyverse)
library(patchwork)
library(ggpubr)

all_contigs_proportion <- read.table(file.choose(), sep = ",", header = TRUE)
#name of the file is chromosome12_contigs.csv
#This contains all the information from Supplemental Table 5, but has been organized to contain only contigs that
#mapped to the sex chromosome of P. picta
#NOTE: Size of Chromosome 12 is 31296.377 kb.

ggscatter(all_contigs_proportion, x = "log2FMcov", y = "proportion",
          color= "segregation", palette = c("lightgrey", "red", "cornflowerblue"), size = 1.5,
          ylab = "Proportion of Chromosome 12", xlab = "log2 F:M Read Depth") + xlim(-2,1.5) +
  geom_vline(xintercept = -0.99,linetype="dotted") + geom_vline(xintercept = 0.6,linetype="dotted")