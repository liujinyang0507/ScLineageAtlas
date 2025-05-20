#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library(ReactomePA)

#library("ggplot2")
#library("DO.db")
#library(R.utils)
#R.utils::setOption("clusterProfiler.download.method",'auto')

pvalueFilter=0.05        
qvalueFilter=0.05         

args<-commandArgs(TRUE)
inputdata1=args[1]

rt=read.table(inputdata1,sep="\t",check.names=F,header=T)
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
write.table(rt,file="id.txt",sep="\t",quote=F,row.names=F)

rt=read.table("id.txt",sep="\t",header=T,check.names=F)   
rt=rt[is.na(rt[,"entrezID"])==F,]

if (ncol(rt)>2){
        geneFC=2^rt$logFC
        gene=rt$entrezID
        names(geneFC)=gene
}else{
        gene=rt$entrezID
}

colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

kk <- enrichPathway(gene, organism = "human", pvalueCutoff = 0.05,
pAdjustMethod = "BH", qvalueCutoff = 0.2, universe, minGSSize = 10,
maxGSSize = 500, readable = T)

GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)

showNum=10
if(nrow(GO)<30){
	showNum=nrow(GO)
}

pdf(file="go_barplot.pdf",width = 9,height = 7)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
		
pdf(file="go_bubble.pdf",width = 9,height = 7)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()
		
if (ncol(rt)>2){
	pdf(file="go_circos.pdf",width = 9,height = 6.5)
	cnet=cnetplot(kk, foldChange=geneFC, showCategory = 5, circular = TRUE, colorEdge = TRUE)
	print(cnet)
	dev.off()
}