rm(list=ls())
library(vegan)
library(ggplot2)
library(ape)

args<-commandArgs(TRUE)
input_file = args[1]
sample = args[2]
out_dir = args[3]

dat <- read.csv(input_file,head = F,sep = " ")
dat <- t(dat)

#计算样本间距离
dis_bray <- dist(dat,method="euclidean")
dis_bray

m="average"
my_hclust <- hclust(dis_bray, method = m)
my_tree <- as.phylo(my_hclust)
labels=paste("clone",seq(ncol(dat)))
my_tree$tip.label <- labels
#write.tree(my_tree, file=paste(out_dir,paste(sample,paste(m,"hclust_lineage.newick",sep = "_"),sep = "_"),sep="/"))
png(file=paste(out_dir,paste(sample,paste(m,"hclust_lineage.png",sep = "_"),sep = "_"),sep="/"))
plot(my_tree)
dev.off()






#层次聚类
# methods <- c("average","centroid","median","complete","single","ward.D2","mcquitty")
# for (m in methods){
#   my_hclust <- hclust(dis_bray, method = m)
#   my_tree <- as.phylo(my_hclust)
#   labels=paste("clone",seq(ncol(dat)))
#   my_tree$tip.label <- labels
#   #write.tree(my_tree, file=paste(out_dir,paste(sample,paste(m,"hclust_lineage.newick",sep = "_"),sep = "_"),sep="/"))
#   png(file=paste(out_dir,paste(sample,paste(m,"hclust_lineage.png",sep = "_"),sep = "_"),sep="/"))
#   plot(my_tree, hang = -1, main = '')
#   dev.off()
# }