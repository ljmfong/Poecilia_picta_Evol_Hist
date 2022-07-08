##### Distribution of Contig Lengths

lengths <- read.table(file.choose(), sep = ",", header = TRUE)
#This is the .csv called contig_lengths.csv
colnames(lengths)[1] <- "scaffold"
dim(lengths)
#[1] 1946   12

lengths$length_kb <- (lengths$length/1000)
lengths$length_mb <- (lengths$length/1000000)
lengths$homo_reads_per_kb <- (lengths$homogametic_median)/(lengths$length_kb)
lengths$hetero_reads_per_kb <- (lengths$heterogametic_median)/(lengths$length_kb)

densityplot <- density(lengths$length_mb)
plot(densityplot, xlab("Contig Lengths (Mb)"), ylab("Density"), main=NULL, xlim=c(0,14))


xlengths_df <- subset(lengths, select = c("scaffold", "segregation", "length", "homogametic_median", "heterogametic_median",
                                          "median_FM_ratio", "length_kb", "homo_reads_per_kb", "hetero_reads_per_kb"))


autosomal_length <- subset(lengths_df, segregation == "autosomal")
xcontig_length <- subset(lengths_df, segregation == "x_contig")
ycontig_length <- subset(lengths_df, segregation == "y_contig")
nodepth_length <-  subset(lengths_df, segregation == "zero")

par(mfrow=c(2,2))

ad <- hist(autosomal_length$homo_reads_per_kb, breaks = 20, col=rgb(0,0,0,1/4),
           xlab = "Homogametic Median Reads per kbp", main = NULL)
ax <- hist(xcontig_length$homo_reads_per_kb, breaks = 20, col=rgb(1,0,0,1/4),
           xlab = "Homogametic Median Reads per kbp", main = NULL)
ay <- hist(ycontig_length$homo_reads_per_kb, breaks = 20, col=rgb(0,0,1,2/4),
           xlab = "Homogametic Median Reads per kbp", main = NULL)
az <- hist(nodepth_length$homo_reads_per_kb, breaks = 20)

plot(ad, col=rgb(0,0,0,1/4), xlab = "Homogametic Median Reads per kbp", main = NULL)
plot(ax, col=rgb(1,0,0,1/4), add = T)
plot(ay, col=rgb(0,0,1,2/4), add = T)
plot(az, col=rgb(0,0,0,1/4), add = T)


md <- hist(autosomal_length$hetero_reads_per_kb, breaks = 20, col=rgb(0,0,0,1/4),
           xlab = "Heterogametic Median Reads per kbp", main = NULL)
mx <- hist(xcontig_length$hetero_reads_per_kb, breaks = 20, col=rgb(1,0,0,1/4),
           xlab = "Heterogametic Median Reads per kbp", main = NULL)
my <- hist(ycontig_length$hetero_reads_per_kb, breaks = 20, col=rgb(0,0,1,2/4),
           xlab = "Heterogametic Median Reads per kbp", main = NULL)
mz <- hist(nodepth_length$hetero_reads_per_kb, breaks = 20,
           xlab = "Heterogametic Median Reads per kbp", main = NULL)

plot(md, col=rgb(0,0,0,1/4), xlab = "Heterogametic Median Reads per kbp", main = NULL)
plot(mx, col=rgb(1,0,0,1/4), add = T)
plot(my, col=rgb(0,0,1,2/4), add = T)
plot(mz, col=rgb(0,0,0,1/4), add = T)


la <- hist(autosomal_length$length_kb, col=rgb(0,0,0,1/4), breaks = 40,
           xlab = "Contig length (kbp)", main = NULL,
           ylim = c(0,1000))
lx <- hist(xcontig_length$length_kb, breaks = 50, col=rgb(1,0,0,1/4),
           xlab = "Contig length (kbp)", main = NULL,
           ylim = c(0,80))
ly <- hist(ycontig_length$length_kb, breaks = 40, col=rgb(0,0,1,2/4),
           xlab = "Contig length (kbp)", main = NULL,
           xlim = c(0,200), ylim = c(0,30))

plot(la, col=rgb(0,0,0,1/4), xlab = "Contig length (kbp)", main = NULL,
     ylim = c(0,1000))
plot(lx, col=rgb(1,0,0,1/4), xlab = "Contig length (kbp)", main = NULL,
     ylim = c(0,80))
plot(ly, col=rgb(0,0,1,2/4), add = T)

lad <- density(autosomal_length$length_kb)
lxd <- density(xcontig_length$length_kb)
lyd <- density(ycontig_length$length_kb)


plot(lad)
plot(lxd)
plot(lyd, xlab = "Y-Contig Length (kbp)",
     col = "blue", main = "Density Plot")

median(ycontig_length$length_kb)
# [1] 22.9295
max(ycontig_length$length_kb)
# [1] 187.103
min(ycontig_length$length_kb)
# [1] 0.896
median(xcontig_length$length_kb)
# [1] 61.176
max(xcontig_length$length_kb)
# [1] 1195.029
min(xcontig_length$length_kb)
# [1] 0.578
median(autosomal_length$length_kb)
# [1] 71.364
max(autosomal_length$length_kb)
# [1] 13339.73
min(autosomal_length$length_kb)
# [1] 0.528


##### Contig Distribution, <10kb removed

rmvd_10kb <- read.table(file.choose(), sep = ",", header = TRUE)
#This is the .csv called rm_10kb.csv
colnames(rmvd_10kb )[1] <- "scaffold"
dim(rmvd_10kb)
#[1] 1276   13
nozeros <- subset(rmvd_10kb, segregation != "zero")
dim(nozeros)
#[1] 1128   13
nozeros$length_mb <- nozeros$length/1000000

densityplot2 <- density(nozeros$length_mb)
plot(densityplot2, xlab("Contig Lengths (Mb)"), ylab("Density"), main=NULL)


rmvd_10kb_df <- subset(rmvd_10kb, select = c("scaffold", "segregation", "homogametic_median", "heterogametic_median",
                                             "median_FM_ratio", "length_kb"))

autosomal_10kblength <- subset(rmvd_10kb_df, segregation == "autosomal")
xcontig_10kblength <- subset(rmvd_10kb_df, segregation == "x_contig")
ycontig_10kblength <- subset(rmvd_10kb_df, segregation == "y_contig")
nodepth_10kblength <-  subset(rmvd_10kb_df, segregation == "zero")


la10kb <- hist(autosomal_10kblength$length_kb, breaks = 50, col=rgb(0,0,0,1/4), 
               xlab = "Contig length (kbp)", main = NULL,
               ylim = c(0,600))
lx10kb <- hist(xcontig_10kblength$length_kb, breaks = 50, col=rgb(1,0,0,1/4),
               xlab = "Contig length (kbp)", main = NULL,
               ylim = c(0,30))
ly10kb <- hist(ycontig_10kblength$length_kb, breaks = 50, col=rgb(0,0,1,2/4),
               xlab = "Contig length (kbp)", main = NULL, 
               xlim = c(0,200), ylim = c(0,30))



plot(lx10kb, col=rgb(1,0,0,1/4), xlab = "Contig length (kbp)", main = NULL,
     ylim = c(0,30))
plot(ly10kb, col=rgb(0,0,1,2/4), add = T)


lyd10kb <- density(ycontig_10kblength$length_kb)
plot(lyd10kb, main = "Density Plot", xlab = "Y-Contig Length (kbp)",
     col = "blue")

lxd10kb <- density(xcontig_10kblength$length_kb)
plot(lxd10kb, main = "Density Plot", xlab = "X-Contig Length (kbp)",
     col = "red", add = T)

lad10kb <- density(autosomal_10kblength$length_kb)
plot(lad10kb, main = "Density Plot", xlab = "Autosomal Contig Length (kbp)")

median(ycontig_10kblength$length_kb)
#[1] 25.757
median(xcontig_10kblength$length_kb)
#[1] 93.95
median(autosomal_10kblength$length_kb)
#[1] 102.7185
