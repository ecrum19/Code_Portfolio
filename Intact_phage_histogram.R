library(ggplot2)
x<-read.table("/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/bact_intact_phage_counts.fna", sep = "\t")
x
my_y <- expression(paste("Number of ", italic("E. coli "), " genomes"))

hist <- ggplot(data=x, aes(V2)) + geom_histogram(bins = 11, fill="grey", col="black") + labs(x="Number of Intact Prophages Identified", y=my_y) + theme_classic() 
hist

install.packages('svglite')
library(svglite)
ggsave(file="/Users/eliascrum/Desktop/prophage_distribution_histogram.svg", plot=hist, width=5, height=5)
