rm(list = ls())
library("treeio")
library("ggtree")
library("ggplot2")
set.seed(123)

args<-commandArgs(TRUE)
input_tree = args[1]
out_dir = args[2]
sample = args[3]

tree_=read.table(input_tree)
tree_=paste(as.character(tree_[1,1]),";")
tree = read.newick(text=tree_)
data = fortify(tree)
dd = subset(data, isTip)
dd = dd[order(dd$y),]
dd$new_label = paste("Clone",c(seq(0,length(dd$label)-1)),sep = "")
for (i in 1:nrow(dd)){
	  tree_=gsub(dd[i,"label"], dd[i,"new_label"], tree_)
} 
tree = read.newick(text=tree_)
write.csv(dd,file=paste(out_dir,paste(sample,"lineage.csv",sep="_"),sep = "/"))
write.tree(tree,file=paste(out_dir,paste(sample,"lineage.newick",sep="_"),sep = "/"))
pdf(file=paste(out_dir,paste(sample,"lineage.png",sep="_"),sep = "/"),onefile = FALSE,
        width = 8,             
	    height =5) 
p = ggtree(tree, layout="rectangular", size=1, col="deepskyblue4",branch.length = 'none') +
	geom_tiplab(hjust = 3,size=5,fontface="italic",offset = 10) +
	geom_point2(aes(subset=!isTip),fill="red",shape=21,size=4) +
	geom_tippoint(size=4, color="deepskyblue4") 
    dev.off()
ggsave(file=paste(out_dir,paste(sample,"lineage.pdf",sep="_"),sep = "/"),p,width = 8, height =6 , dpi = 300)
