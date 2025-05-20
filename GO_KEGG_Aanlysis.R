library(clusterProfiler)
library(ggplot2)
library("org.Hs.eg.db")
library(forcats)
library(GO.db, lib.loc="/home/ubuntu/miniconda3/envs/CellLineage2/lib/R/library")
library(dplyr)
library(stringr)
options(clusterProfiler.download.method = "wget")
pvalueFilter=0.05         #p值过滤条件
qvalueFilter=0.05         #矫正后的p值过滤条件

args<-commandArgs(TRUE)
input_file = args[1]
out_dir = args[2]
sample = args[3]
df_sig <- read.csv(input_file, header = T)
group <- data.frame(gene=df_sig$gene,
		    group=df_sig$label)

Gene_ID <- bitr(df_sig$gene, fromType="SYMBOL",
		toType="ENTREZID",OrgDb="org.Hs.eg.db")
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
data_GO <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

data_GO_sim <- simplify(data_GO, 
                        cutoff=0.7, 
                        by="p.adjust", 
                        select_fun=min)

data_GO_sim_fil <- data_GO_sim@compareClusterResult
data_GO_sim_fil = data_GO_sim_fil[is.na(data_GO_sim_fil[,"pvalue"])==F,] 
data_GO_sim_fil = data_GO_sim_fil[(data_GO_sim_fil$pvalue<pvalueFilter & data_GO_sim_fil$qvalue<qvalueFilter),]
write.csv(data_GO_sim_fil,paste(out_dir,paste(sample,"GO_total.csv",sep="_"),sep = "/"),row.names=FALSE,quote=FALSE)

data_GO_sim_fil_plot = data_GO_sim_fil %>% group_by(group,ONTOLOGY) %>% filter(row_number() <= 5)
data_GO_sim_fil_plot  = data_GO_sim_fil_plot [is.na(data_GO_sim_fil_plot[,"qvalue"])==F,]

data_GO_sim_fil_plot$term_ <- paste(data_GO_sim_fil_plot$ID, data_GO_sim_fil_plot$ONTOLOGY, sep = ': ')
data_GO_sim_fil_plot$term <- paste(data_GO_sim_fil_plot$term_, data_GO_sim_fil_plot$Description, sep = ': ')
data_GO_sim_fil_plot <- select(data_GO_sim_fil_plot,-term_)
write.csv(data_GO_sim_fil_plot,paste(out_dir,paste(sample,"GO_visualization.csv",sep="_"),sep = "/"),row.names=FALSE,quote=FALSE)
head(data_GO_sim_fil_plot)
P <- ggplot(data_GO_sim_fil_plot,
             aes(y=term,x=group))+
  geom_point(aes(size=Count,color=pvalue))+   
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Group",y="GO term") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=50)) +
  theme(plot.margin=unit(rep(1,4),'cm')) +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 12)) +
  theme_bw()
ggsave(paste(out_dir,paste(sample,"GO.png",sep="_"),sep = "/"), P, width = 6, height = 10, dpi = 300)
ggsave(paste(out_dir,paste(sample,"GO.pdf",sep="_"),sep = "/"), P, width = 6, height = 10, dpi = 300)

data_KEGG <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichKEGG",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

data_KEGG_sim_fil <- data_KEGG@compareClusterResult
data_KEGG_sim_fil = data_KEGG_sim_fil[is.na(data_KEGG_sim_fil[,"pvalue"])==F,] 
data_KEGG_sim_fil = data_KEGG_sim_fil[(data_KEGG_sim_fil$pvalue<pvalueFilter & data_KEGG_sim_fil$qvalue<qvalueFilter),]
write.csv(data_KEGG_sim_fil,paste(out_dir,paste(sample,"KEGG_total.csv",sep="_"),sep = "/"),row.names=FALSE,quote=FALSE)

data_KEGG_sim_fil_plot = data_KEGG_sim_fil %>% group_by(group) %>% filter(row_number() <= 5)
data_KEGG_sim_fil_plot  = data_KEGG_sim_fil_plot [is.na(data_KEGG_sim_fil_plot[,"qvalue"])==F,]
write.csv(data_KEGG_sim_fil_plot,paste(out_dir,paste(sample,"KEGG_visualization.csv",sep="_"),sep = "/"),row.names=FALSE,quote=FALSE)

#data_KEGG_sim_fil_plot=read.csv(paste(out_dir,"kegg_visualization.csv",sep = "/"))
data_KEGG_sim_fil_plot$term <- paste(data_KEGG_sim_fil_plot$ID, data_KEGG_sim_fil_plot$Description, sep = ': ')
p = ggplot(data_KEGG_sim_fil_plot,
            aes(y=term,x=group))+
  geom_point(aes(size=Count,color=pvalue))+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Group",y="KEGG term") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=30)) +
  theme(plot.margin=unit(rep(1,4),'cm')) +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 12))+
  theme_bw()
ggsave(paste(out_dir,paste(sample,"KEGG.png",sep="_"),sep="/"), p, width = 6, height = 10, dpi = 300)
ggsave(paste(out_dir,paste(sample,"KEGG.pdf",sep="_"),sep="/"), p, width = 6, height = 10, dpi = 300)