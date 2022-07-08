#Kmers

library(ggplot2)

#### Plotting only the histo counts after histo command from jellyfish of pooled samples.
#Sexes pooled separately
#No change in any of the files (files have been trimmed already)
#This is normalized to a coverage of 30x

male_kmer_counts <- read.table(file.choose(), header=FALSE, sep=" ")
#file is all_male_DNA_kmer_norm30.histo

plot(male_kmer_counts, col="darkorchid3", type='l', ylim=c(0,1e7), xlim=c(0,120), 
     xlab="k-mer Coverage", ylab="No. of k-mers counted",
     cex.lab=0.7, cex.axis = 0.7, xaxt = "n")
axis(1, at = seq(0,120, by = 20), cex.axis=0.7)

#Now read in the females
female_kmer_counts <- read.table(file.choose(), header=FALSE, sep=" ")
#file is all_female_DNA_kmer_norm30.histo


lines(female_kmer_counts, type= 'l', col="red")


#### Plotting only the histo counts after histo command from jellyfish of pooled samples.

















###################
counts<-read.csv(file.choose(),header=TRUE)
#file was named picta_kmer_counts.csv
colnames(counts)[1] <- "coverage"

first_100 <- subset(counts, coverage>0 & coverage<101)

ggplot(first_100, aes(x=coverage,y=count,fill=sex)) +
  geom_bar(stat = "identity", alpha=0.5) +
  ylab("No. of k-mers counted") + ylim(0,6e6) +
  xlab("k-mer coverage")



first_80 <- subset(counts, coverage>0 & coverage<81)

ggplot(first_80, aes(x=coverage,y=count,fill=sex)) +
  geom_bar(stat = "identity", alpha=0.5) +
  ylab("No. of k-mers counted") + ylim(0,6e6) +
  xlab("k-mer coverage")

g3_l80<- subset(counts, coverage>3 & coverage<81)

ggplot(g3_l80, aes(x=coverage,y=count,fill=sex)) +
  geom_bar(stat = "identity", alpha=0.5) +
  ylab("No. of k-mers counted") + ylim(0,6e6) +
  xlab("k-mer coverage")



########### OLD ############

bins<-first_4000$bin
only_fem<-first_4000$log_fem
only_mal<-first_4000$log_mal

fem_set <- data.frame(cbind(bins, only_fem))
mal_set <- data.frame(cbind(bins, only_mal))

ggplot(fem_set, aes(x=bins, y=only_fem)) + 
  geom_bar(stat = "identity") + ylab("Average Female Counts (Log)")

ggplot(mal_set, aes(x=bins, y=only_mal)) + 
  geom_bar(stat = "identity") + ylab("Average Male Counts (Log)")

first_500<-counts[1:500,]

bins<-first_500$bin
only_fem<-first_500$log_fem
only_mal<-first_500$log_mal

fem_set <- data.frame(cbind(bins, only_fem))
mal_set <- data.frame(cbind(bins, only_mal))

ggplot(fem_set, aes(x=bins, y=only_fem)) + 
  geom_bar(stat = "identity") + ylab("Average Female Counts (Log)") +
  xlab("k-mer coverage")


ggplot(mal_set, aes(x=bins, y=only_mal)) + 
  geom_bar(stat = "identity") + ylab("Average Male Counts (Log)") +
  xlab("k-mer coverage")

first_300 <- counts[1:300,]

bins <- first_300$bin
only_fem <-first_300$log_fem
only_mal <-first_300$log_mal

fem_set <- data.frame(cbind(bins, only_fem))
mal_set <- data.frame(cbind(bins, only_mal))

#Base plot comparison

plot(density(fem_set$only_fem))
lines(density(mal_set$only_mal))

#Pretty plot overlay

combo_set <- data.frame(cbind(fem_set,mal_set$only_mal))
colnames(combo_set)[3] <- "only_mal"
max(combo_set$only_fem)
#[1] 7.459474
max(combo_set$only_mal)
#[1] 7.46306

ggplot() + 
  geom_area(data = combo_set, aes(x = as.numeric(row.names(combo_set)),y=only_fem),
            fill = 'darkorange', alpha = 0.4) +
  geom_area(data = combo_set, aes(x=as.numeric(row.names(combo_set)),y=only_mal),
            fill = 'grey94') +
  geom_line(data = fem_set, aes(x=as.numeric(row.names(combo_set)),y=only_fem), 
               stat = 'identity', colour = 'darkorange') +
  geom_line(data = mal_set, aes(x=as.numeric(row.names(combo_set)),y=only_mal), 
               stat = 'identity', colour = 'darkorchid') +
  theme_classic() + xlab("k-mer coverage") + ylab("k-mers counted")



#### Write the script for the actual data ####

library(ggplot2)

#Read the file in - it's a big file, headers will include
#Kmer MaleCoverage  FemaleCoverage

#data file:

combo_set<-read.table()

ggplot() + 
  geom_area(data = combo_set, aes(x=as.numeric(row.names(combo_set)),y=FemaleCoverage),
            fill = 'darkorange', alpha = 0.4) +
  geom_area(data = combo_set, aes(x=as.numeric(row.names(combo_set)),y=MaleCoverage),
            fill = 'grey94') +
  geom_line(data = fem_set, aes(x=as.numeric(row.names(combo_set)),y=FemaleCoverage), 
            stat = 'identity', colour = 'darkorange') +
  geom_line(data = mal_set, aes(x=as.numeric(row.names(combo_set)),y=MaleCoverage), 
            stat = 'identity', colour = 'darkorchid') +
  theme_classic() + xlab("k-mer coverage") + ylab("k-mers counted")

