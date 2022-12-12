#Male-Biased Gene Expression Script for P. picta

library(DESeq)
library(DESeq2)
library(edgeR)
library(devtools)
library(EnhancedVolcano)
library(limma)
library(DEFormats)
library(ggplot2)

#input expression data - taken from salmon counts

raw_gene_counts <- read.table(file.choose(), header = TRUE, sep = ",") #This was my salmon_counts.csv file
colnames(raw_gene_counts)[1] <- "gene_id"
gene_names_row <- raw_gene_counts[,-1]
rownames(gene_names_row) <- raw_gene_counts[,1]

gene_expression <- data.frame(gene_names_row) #this is now a matrix
#From the visualizing gene expression Rscript

conditions<-factor(c("F","F","F","F","F","F","M","M","M","M","M","M",
                     "F","F","F","F","F","F","M","M","M","M","M","M",
                     "F","F","F","F","F","F","M","M","M","M","M","M",
                     "F","F","F","F","F","F","M","M","M","M","M","M"))
length(conditions)
dim(gene_expression)
#You will want your conditions to match the number of labelled columns 
#i.e. the length should match the length of your dataframe/matrix

all_gene_expression <- DGEList(counts = gene_expression,group=conditions)
plotMDS(all_gene_expression)

plotSmear(all_gene_expression,pair=NULL, de.tags=NULL, xlab="Average logCPM", ylab="logFC", pch=19, cex=0.2, smearWidth=0.5, panel.first=grid(), smooth.scatter=FALSE, lowess=FALSE)
cpm_all_gene_expression_raw <- cpm(all_gene_expression)

sample1 <- density(log2(cpm_all_gene_expression_raw[,1]))
sample2 <- density(log2(cpm_all_gene_expression_raw[,2]))
sample3 <- density(log2(cpm_all_gene_expression_raw[,3]))
sample4 <- density(log2(cpm_all_gene_expression_raw[,4]))
sample5 <- density(log2(cpm_all_gene_expression_raw[,5]))
sample6 <- density(log2(cpm_all_gene_expression_raw[,6]))
sample7 <- density(log2(cpm_all_gene_expression_raw[,7]))
sample8 <- density(log2(cpm_all_gene_expression_raw[,8]))
sample9 <- density(log2(cpm_all_gene_expression_raw[,9]))
sample10 <- density(log2(cpm_all_gene_expression_raw[,10]))
sample11 <- density(log2(cpm_all_gene_expression_raw[,11]))
sample12 <- density(log2(cpm_all_gene_expression_raw[,12]))
sample13 <- density(log2(cpm_all_gene_expression_raw[,13]))
sample14 <- density(log2(cpm_all_gene_expression_raw[,14]))
sample15 <- density(log2(cpm_all_gene_expression_raw[,15]))
sample16 <- density(log2(cpm_all_gene_expression_raw[,16]))
sample17 <- density(log2(cpm_all_gene_expression_raw[,17]))
sample18 <- density(log2(cpm_all_gene_expression_raw[,18]))
sample19 <- density(log2(cpm_all_gene_expression_raw[,19]))
sample20 <- density(log2(cpm_all_gene_expression_raw[,20]))
sample21 <- density(log2(cpm_all_gene_expression_raw[,21]))
sample22 <- density(log2(cpm_all_gene_expression_raw[,22]))
sample23 <- density(log2(cpm_all_gene_expression_raw[,23]))
sample24 <- density(log2(cpm_all_gene_expression_raw[,24]))
sample25 <- density(log2(cpm_all_gene_expression_raw[,25]))
sample26 <- density(log2(cpm_all_gene_expression_raw[,26]))
sample27 <- density(log2(cpm_all_gene_expression_raw[,27]))
sample28 <- density(log2(cpm_all_gene_expression_raw[,28]))
sample29 <- density(log2(cpm_all_gene_expression_raw[,29]))
sample30 <- density(log2(cpm_all_gene_expression_raw[,30]))
sample31 <- density(log2(cpm_all_gene_expression_raw[,31]))
sample32 <- density(log2(cpm_all_gene_expression_raw[,32]))
sample33 <- density(log2(cpm_all_gene_expression_raw[,33]))
sample34 <- density(log2(cpm_all_gene_expression_raw[,34]))
sample35 <- density(log2(cpm_all_gene_expression_raw[,35]))
sample36 <- density(log2(cpm_all_gene_expression_raw[,36]))
sample37 <- density(log2(cpm_all_gene_expression_raw[,37]))
sample38 <- density(log2(cpm_all_gene_expression_raw[,38]))
sample39 <- density(log2(cpm_all_gene_expression_raw[,39]))
sample40 <- density(log2(cpm_all_gene_expression_raw[,40]))
sample41 <- density(log2(cpm_all_gene_expression_raw[,41]))
sample42 <- density(log2(cpm_all_gene_expression_raw[,42]))
sample43 <- density(log2(cpm_all_gene_expression_raw[,43]))
sample44 <- density(log2(cpm_all_gene_expression_raw[,44]))
sample45 <- density(log2(cpm_all_gene_expression_raw[,45]))
sample46 <- density(log2(cpm_all_gene_expression_raw[,46]))
sample47 <- density(log2(cpm_all_gene_expression_raw[,47]))
sample48 <- density(log2(cpm_all_gene_expression_raw[,48]))
plot(sample1, xlab = "CPM", ylab = "Density",type="l",lwd=2,main="Raw log2 cpm")
lines(sample2, type="l",lwd=2)
lines(sample3, type="l",lwd=2)
lines(sample4, type="l",lwd=2)
lines(sample5, type="l",lwd=2)
lines(sample6, type="l",lwd=2)
lines(sample7, type="l",lwd=2)
lines(sample8, type="l",lwd=2)
lines(sample9, type="l",lwd=2)
lines(sample10, type="l",lwd=2)
lines(sample11, type="l",lwd=2)
lines(sample12, type="l",lwd=2)
lines(sample13, type="l",lwd=2)
lines(sample14, type="l",lwd=2)
lines(sample15, type="l",lwd=2)
lines(sample16, type="l",lwd=2)
lines(sample17, type="l",lwd=2)
lines(sample18, type="l",lwd=2)
lines(sample19, type="l",lwd=2)
lines(sample20, type="l",lwd=2)
lines(sample21, type="l",lwd=2)
lines(sample22, type="l",lwd=2)
lines(sample23, type="l",lwd=2)
lines(sample24, type="l",lwd=2)
lines(sample25, type="l",lwd=2)
lines(sample26, type="l",lwd=2)
lines(sample27, type="l",lwd=2)
lines(sample28, type="l",lwd=2)
lines(sample29, type="l",lwd=2)
lines(sample30, type="l",lwd=2)
lines(sample31, type="l",lwd=2)
lines(sample32, type="l",lwd=2)
lines(sample33, type="l",lwd=2)
lines(sample34, type="l",lwd=2)
lines(sample35, type="l",lwd=2)
lines(sample36, type="l",lwd=2)
lines(sample37, type="l",lwd=2)
lines(sample38, type="l",lwd=2)
lines(sample39, type="l",lwd=2)
lines(sample40, type="l",lwd=2)
lines(sample41, type="l",lwd=2)
lines(sample42, type="l",lwd=2)
lines(sample43, type="l",lwd=2)
lines(sample44, type="l",lwd=2)
lines(sample45, type="l",lwd=2)
lines(sample46, type="l",lwd=2)
lines(sample47, type="l",lwd=2)
lines(sample48, type="l",lwd=2)
#Check normalization of count data

expr_genes <- DGEList(counts=all_gene_expression,group=conditions)
genes_norm <- calcNormFactors(expr_genes)
plotMDS(genes_norm)
plotSmear(genes_norm,pair=NULL, de.tags=NULL, xlab="Average logCPM", ylab="logFC", pch=19, cex=0.2, smearWidth=0.5, panel.first=grid(), smooth.scatter=FALSE, lowess=FALSE)
cpm_genes_norm <- cpm(genes_norm)
sample1 <- density(log2(cpm_genes_norm[,1]))
sample2 <- density(log2(cpm_genes_norm[,2]))
sample3 <- density(log2(cpm_genes_norm[,3]))
sample4 <- density(log2(cpm_genes_norm[,4]))
sample5 <- density(log2(cpm_genes_norm[,5]))
sample6 <- density(log2(cpm_genes_norm[,6]))
sample7 <- density(log2(cpm_genes_norm[,7]))
sample8 <- density(log2(cpm_genes_norm[,8]))
sample9 <- density(log2(cpm_genes_norm[,9]))
sample10 <- density(log2(cpm_genes_norm[,10]))
sample11 <- density(log2(cpm_genes_norm[,11]))
sample12 <- density(log2(cpm_genes_norm[,12]))
sample13 <- density(log2(cpm_genes_norm[,13]))
sample14 <- density(log2(cpm_genes_norm[,14]))
sample15 <- density(log2(cpm_genes_norm[,15]))
sample16 <- density(log2(cpm_genes_norm[,16]))
sample17 <- density(log2(cpm_genes_norm[,17]))
sample18 <- density(log2(cpm_genes_norm[,18]))
sample19 <- density(log2(cpm_genes_norm[,19]))
sample20 <- density(log2(cpm_genes_norm[,20]))
sample21 <- density(log2(cpm_genes_norm[,21]))
sample22 <- density(log2(cpm_genes_norm[,22]))
sample23 <- density(log2(cpm_genes_norm[,23]))
sample24 <- density(log2(cpm_genes_norm[,24]))
sample1 <- density(log2(cpm_genes_norm[,25]))
sample2 <- density(log2(cpm_genes_norm[,26]))
sample3 <- density(log2(cpm_genes_norm[,27]))
sample4 <- density(log2(cpm_genes_norm[,28]))
sample5 <- density(log2(cpm_genes_norm[,29]))
sample6 <- density(log2(cpm_genes_norm[,30]))
sample7 <- density(log2(cpm_genes_norm[,31]))
sample8 <- density(log2(cpm_genes_norm[,32]))
sample9 <- density(log2(cpm_genes_norm[,33]))
sample10 <- density(log2(cpm_genes_norm[,34]))
sample11 <- density(log2(cpm_genes_norm[,35]))
sample12 <- density(log2(cpm_genes_norm[,36]))
sample13 <- density(log2(cpm_genes_norm[,37]))
sample14 <- density(log2(cpm_genes_norm[,38]))
sample15 <- density(log2(cpm_genes_norm[,39]))
sample16 <- density(log2(cpm_genes_norm[,40]))
sample17 <- density(log2(cpm_genes_norm[,41]))
sample18 <- density(log2(cpm_genes_norm[,42]))
sample19 <- density(log2(cpm_genes_norm[,43]))
sample20 <- density(log2(cpm_genes_norm[,44]))
sample21 <- density(log2(cpm_genes_norm[,45]))
sample22 <- density(log2(cpm_genes_norm[,46]))
sample23 <- density(log2(cpm_genes_norm[,47]))
sample24 <- density(log2(cpm_genes_norm[,48]))
plot(sample1, xlab = "CPM", ylab = "Density",type="l",lwd=2,main="Norm log2 cpm")
lines(sample2, type="l",lwd=2)
lines(sample3, type="l",lwd=2)
lines(sample4, type="l",lwd=2)
lines(sample5, type="l",lwd=2)
lines(sample6, type="l",lwd=2)
lines(sample7, type="l",lwd=2)
lines(sample8, type="l",lwd=2)
lines(sample9, type="l",lwd=2)
lines(sample10, type="l",lwd=2)
lines(sample11, type="l",lwd=2)
lines(sample12, type="l",lwd=2)
lines(sample13, type="l",lwd=2)
lines(sample14, type="l",lwd=2)
lines(sample15, type="l",lwd=2)
lines(sample16, type="l",lwd=2)
lines(sample17, type="l",lwd=2)
lines(sample18, type="l",lwd=2)
lines(sample19, type="l",lwd=2)
lines(sample20, type="l",lwd=2)
lines(sample21, type="l",lwd=2)
lines(sample22, type="l",lwd=2)
lines(sample23, type="l",lwd=2)
lines(sample24, type="l",lwd=2)
lines(sample25, type="l",lwd=2)
lines(sample26, type="l",lwd=2)
lines(sample27, type="l",lwd=2)
lines(sample28, type="l",lwd=2)
lines(sample29, type="l",lwd=2)
lines(sample30, type="l",lwd=2)
lines(sample31, type="l",lwd=2)
lines(sample32, type="l",lwd=2)
lines(sample33, type="l",lwd=2)
lines(sample34, type="l",lwd=2)
lines(sample35, type="l",lwd=2)
lines(sample36, type="l",lwd=2)
lines(sample37, type="l",lwd=2)
lines(sample38, type="l",lwd=2)
lines(sample39, type="l",lwd=2)
lines(sample40, type="l",lwd=2)
lines(sample41, type="l",lwd=2)
lines(sample42, type="l",lwd=2)
lines(sample43, type="l",lwd=2)
lines(sample44, type="l",lwd=2)
lines(sample45, type="l",lwd=2)
lines(sample46, type="l",lwd=2)
lines(sample47, type="l",lwd=2)
lines(sample48, type="l",lwd=2)

#Normalise and extract rpkm
#Gene length is taken from the gene_length found in RSEM.isoform.results file
expr_genes <- DGEList(counts=all_gene_expression)
genes_norm <- calcNormFactors(expr_genes)
gene_length_all <- read.table(file.choose(), stringsAsFactors = FALSE) #this should be a file with two columns 
#one with all the gene names, and the other with the gene lengths: there is no header (just V1 and V2)
#Always check your gene_length_all file in case something looks wrong
#I saved it as a tab-delimited file and this helped
#The file is called salmon_gene_length.txt
expressed_genes<-rownames(all_gene_expression)
length(expressed_genes)
gene_length <- subset(gene_length_all, V1 %in% expressed_genes)
gene_length_vector <- c(gene_length$V2)
all(gene_length$V1 == rownames(all_gene_expression)) #should print TRUE; we want the row names to match


rpkm_norm <- rpkm(genes_norm, log=FALSE, gene.length=gene_length_vector)
write.table(rpkm_norm, file="all_gene_expression_rpkm.txt",quote=F, sep="\t")
plotMD(rpkm_norm)
dim(rpkm_norm)


library(pheatmap)
library(pvclust)
library(colorRamps)
palette2 <-colorRamps::"matlab.like2"(n=200)
bootstraps<- pvclust(log2(rpkm_norm+1), method.hclust="average", method.dist="euclidean")
plot(bootstraps)

pheatmap(log2(rpkm_norm+1), border_color=NA,show_colnames=T,show_rownames=T,color = colorRampPalette(c("gray30","gray48","gray68","white","yellow","yellow2","gold2"))(100),clustering_distance_cols = "euclidean", clustering_method="average")

#We have a total of 25,378 genes - we now want to find:
#1) Expressed only in females
#2) Expressed only in males
#3) Ratio of those genes present in both
#4) Genes present only in males not in females
### Filtering the genes:
#When needing to pull from the file that you had already saved and not wanting to redo all the steps previously:

rpkm_norm_df <- read.table(file.choose(), header = TRUE)

rpkm_norm_df <- data.frame(rpkm_norm)

female_genes <- rpkm_norm_df[,c("T1_Dame","T1_F1","T1_F2","T1_F3", "T1_F5", "T1_F6",
                                "T13_Dame", "T13_F1","T13_F3","T13_F4","T13_F6","T13_F7",
                                "G3_Dame", "G3_F1", "G3_F2", "G3_F3", "G3_F4", "G3_F5",
                                "G5_Dame", "G5_F1","G5_F2", "G5_F3", "G5_F4", "G5_F6")]
female_genes$female_median <- apply(female_genes[,-1], 1, median)


male_genes <- rpkm_norm_df[,c("T1_Sire","T1_M3","T1_M5","T1_M6","T1_M7","T1_M8", 
                              "T13_Sire","T13_M2","T13_M3","T13_M5","T13_M7","T13_M9",
                              "G3_Sire", "G3_M3", "G3_M6", "G3_M7", "G3_M10", "G3_M11",
                              "G5_Sire", "G5_M1", "G5_M2", "G5_M4", "G5_M5", "G5_M6")]
male_genes$male_median <- apply(male_genes[,-1], 1, median)

rpkm_medians <- cbind(female_genes,male_genes)
rpkm_medians_df <- rpkm_medians[c(25,50)]
rpkm_medians_df$log2FC <- log2((rpkm_medians_df$female_median+0.01)/(rpkm_medians_df$male_median+0.01))


#Remove genes that have a median of zero (this will filtre out the lowly expressed genes)
library(tibble)
adjust_rows <- tibble::rownames_to_column(rpkm_medians_df, "gene")

library(dplyr)

rpkm_medians_HE <- adjust_rows %>%
  mutate(female_median == 0 & male_median == 0)



rpkm_medians_HE_filtered <- subset(rpkm_medians_HE, `female_median == 0 & male_median == 0` == "FALSE")
dim(rpkm_medians_HE_filtered)
#[1] 22493     5 #A total of 22,493

#1) Expressed only in females
female_uniq_express <- subset(rpkm_medians_HE_filtered, male_median == 0)
dim(female_uniq_express)
#[1] 671   5 # 671 uniquely expressed genes
write.table(female_uniq_express, file="female_uniq_gene_expression_rpkm.txt",quote=F, sep="\t")


#2) Expressed only in males
male_uniq_express <- subset(rpkm_medians_HE_filtered, female_median == 0)
dim(male_uniq_express)
#[1] 353   5 #353 uniquely expressed genes
write.table(male_uniq_express, file="male_uniq_gene_expression_rpkm.txt",quote=F, sep="\t")

#3) Ratio of those genes present in both
write.table(rpkm_medians_HE_filtered, file="HE_gene_expression_rpkm.txt",quote=F, sep="\t")
#The rest was done in excel
#I removed the duplicated genes (e.g. male_ or female_uniq_express)

#Plot
ratio_present_both <- read.csv(file.choose(),header = TRUE, sep =",")
#name of file is DGE_genes_biases.csv
colnames(ratio_present_both)[1] <- "gene"
uniq_male_median <- male_uniq_express[1:4]
uniq_male_median$expression_type <- c("male-bias")
uniq_female_median <- female_uniq_express[1:4]
uniq_female_median$expression_type <- c("female-bias")

appended_plot <- rbind(ratio_present_both,uniq_female_median,uniq_male_median)
appended_plot$label <- ifelse(appended_plot$expression_type == "unbiased", "background", "point")

ggplot(appended_plot) +  
  geom_point(aes(x = log2(male_median+1), y = log2(female_median+1), colour = expression_type)) +
  geom_point(data=subset(appended_plot, label == 'point'),
             (aes(x = log2(male_median+1), y = log2(female_median+1), colour = expression_type))) +
  scale_color_manual(values=c("orange", "purple", "lightgrey")) +
  xlab("Log2 median male counts (RPKM)") + 
  ylab("Log2 median female counts (RPKM)") + theme_classic()


########## Checking my male-biased genes

library(ggplot2)
library(dplyr)

#We can look at the number of genes that are considered male biased (FC > 2)
#The log2FC calculated by R is very different than the one calculated by excel.
#Therefore, I did it on excel and am doing the rest of the analysis from the excel sheet:
checking_male_bias <- read.table(file.choose(), header = TRUE, sep = ",")
#File is: C:\Users\17802\Documents\Lydia\Picta_project\male_limited_expression\Updated_DGE\checking_male_biased_genes.csv

colnames(checking_male_bias)[1] <- "gene"

#This file has all of the genes in RPKM

rpkm_medians_HE <- checking_male_bias %>%
  mutate(female_median_rpkm == 0 & male_median_rpkm == 0)

rpkm_medians_HE_filtered <- subset(rpkm_medians_HE, `female_median_rpkm == 0 & male_median_rpkm == 0` == "FALSE")
dim(rpkm_medians_HE_filtered)
#[1] 22715     6 #A total of 22,715

#1) Expressed only in females
female_uniq_express <- subset(rpkm_medians_HE_filtered, male_median_rpkm == 0)
dim(female_uniq_express)
#[1] 564   6 # 564 uniquely expressed genes
write.table(female_uniq_express, file="female_uniq_gene_expression_rpkm_excel.txt",quote=F, sep="\t")

#1.5) Expressed only in females and HIGHLY expressed (RPKM >2)
female_uniq_express_HE <- subset(female_uniq_express, female_median_rpkm > 2)
dim(female_uniq_express_HE)

#2) Expressed only in males
male_uniq_express <- subset(rpkm_medians_HE_filtered, female_median_rpkm == 0)
dim(male_uniq_express)
#[1] 317   6 #317 uniquely expressed genes
write.table(male_uniq_express, file="male_uniq_gene_expression_rpkm_excel.txt",quote=F, sep="\t")

#2.5) Expressed only in males and HIGHLY expressed (RPKM >2)
male_uniq_express_HE <- subset(male_uniq_express, male_median_rpkm > 2)
dim(male_uniq_express_HE)

#To note: Few of the unique genes are highly expressed (RPKM >2)

#3) Ratio of those genes present in both, but are male biased (log2FC > 2)

male_biased_express <- subset(rpkm_medians_HE_filtered, log2FC_MF > 2)
dim(male_biased_express)
#[1] 383   6 #We have 383 male-biased expressed genes, including those that are male unique

male_biased_express_HE <- subset(male_biased_express, male_median_rpkm >2)
dim(male_biased_express_HE)
#[1] 27  6
# Of which, only 27 of them can be considered highly expressed, despite the low expression in females
write.table(male_biased_express_HE, file="male_biased_expression_excel.txt",quote=F, sep="\t")

#To plot the filtered genes:



ggplot(rpkm_medians_HE_filtered, aes(x = male_median_rpkm, y = female_median_rpkm)) + 
  geom_point() + xlab("median male counts (RPKM)") + ylab("median female counts (RPKM)") +
  geom_point(data = male_biased_express_HE, aes(x = male_median_rpkm, y = female_median_rpkm),
             colour = 'red') + theme_classic()

#Log plot
ggplot(rpkm_medians_HE_filtered, aes(x = log(male_median_rpkm+1), y = log(female_median_rpkm+1))) + 
  geom_point() + xlab("log median male counts (RPKM)") + ylab("log median female counts (RPKM)") +
  geom_point(data = male_biased_express_HE, aes(x = log(male_median_rpkm+1), y = log(female_median_rpkm+1)),
             colour = 'red') + theme_classic()