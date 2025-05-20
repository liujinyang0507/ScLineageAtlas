# install.packages("devtools")
# library(devtools)
# install_github("lhe17/nebula")
# install.packages('randomcoloR')

#rm(list=ls())
library(Seurat)
library(SeuratDisk)
library(nebula)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(showtext)
library(scales)
library(randomcoloR)

args<-commandArgs(TRUE)
h5ad_file <- args[1]
obs_file  <- args[2]
patient_file <- args[3]
clones_dir <- args[4]
out_dir <- args[5]

# Conver .h5ad to .h5seurat
Convert(h5ad_file,dest = 'h5seurat',overwrite = T,assay = 'RNA')

# Load h5seurat file
h5seurat <- gsub("h5ad","h5seurat",h5ad_file)
adata <- LoadH5Seurat(h5seurat,meta.data=F)

# Load meta.data (obs from h5ad)
meta <- read.csv(obs_file,row.names = 1)

# AddMetadata to Seurat Object
adata <- AddMetaData(adata,meta)
print("**********adata infor**********")
head(adata)

# Load Patient annotation
pat <- read.csv(patient_file)
pat <- unique(pat)

# Convert patient annotation to a vector
patient <- as.vector(pat$sample_title)
#names(patient) <- pat$run_accession


plot_multiple_volcano <- function(output_sig,out_dir){
  # define a color palette
  dfcol<-data.frame(x=c(1:length(unique(output_sig$label))),
                    y=0,
                    label = unique(output_sig$label))
  mycol <- distinctColorPalette(length(unique(output_sig$label)))
  
  top10_sig <- data.frame()
  for (i in unique(output_sig$label)) {
    output_sub <- output_sig[output_sig$label==i,]
    top10 <- output_sub %>% top_n(10,abs(logFC))
    top10_sig <- rbind(top10_sig,top10)
  }
  
  # 多组火山图
  p1 <- ggplot()+
    geom_jitter(data = output_sig,
                aes(x = label, y = logFC, color = expression),
                size = 0.6,
                width = 0.3) +
    # scale_y_reverse()+
    scale_y_continuous(limits = c(10,1e18),expand=c(0,0))+
    scale_color_manual(values = c("Down-regulated"='blue',"Up-regulated"='red'))+
    geom_text_repel(
      data=top10_sig,
      aes(x=label,y=logFC,label=gene),
      # force = 1.2,
      # arrow = arrow(length = unit(0.008, "npc"),
      #               type = "open", ends = "last"),
      size=3
    )+
    # geom_col(data = dfbar,
    #          mapping = aes(x = x,y = y),
    #          fill = "#dcdcdc",alpha = 0.6)+
    # geom_col(data = dfbar1,
    #          mapping = aes(x = x,y = y),
    #          fill = "#dcdcdc",alpha = 0.6)+
    theme_minimal()+
    theme(
      axis.title = element_text(size = 13,
                                color = "black",
                                face = "bold"),
      axis.line.y.left = element_line(color = "black",
                                      size = 1.2),
      axis.line.y = element_blank(),
      # axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = c(1,0),
      legend.text = element_text(size = 10),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # axis.text.y = element_blank(),
    )+
    ggtitle(p)
  # coord_flip()
  # ylim(10,1e18)
  p1
  
  p2 <- ggplot()+
    # geom_col(data = dfbar,
    #          mapping = aes(x = x,y = y),
    #          fill = "#dcdcdc",alpha = 0.6)+
    # geom_col(data = dfbar1,
    #          mapping = aes(x = x,y = y),
    #          fill = "#dcdcdc",alpha = 0.6)+
    geom_jitter(data = output_sig,
                aes(x = label, y = logFC, color = expression),
                size = 0.6,
                width =0.3) +
    # ylim(-10,10) +
    geom_tile(data = dfcol,
              aes(x=x,y=y),
              height=3,
              color = "black",
              fill = mycol,
              alpha = 1,
              show.legend = F)+
    # labs(x="Cluster")+
    geom_text(data=dfcol,
              aes(x=x,y=y,label=label),
              size = 5,
              color ="white")+
    scale_y_continuous(limits = c(-10,10),expand=c(0,0))+
    scale_color_manual(values = c("Down-regulated"='blue',"Up-regulated"='red'))+
    geom_text_repel(
      data=top10_sig,
      aes(x=label,y=logFC,label=gene),
      force = 2,
      # max.overlaps = 10,
      size=3,
      # arrow = arrow(length = unit(0.008, "npc"),
      #               type = "open", ends = "last")
    )+
    theme_minimal()+
    theme(
      axis.title = element_text(size = 13,
                                color = "black",
                                face = "bold"),
      axis.line.y.left = element_line(color = "black",
                                      size = 1.2),
      axis.line.y = element_blank(),
      # axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = c(1,0),
      legend.text = element_text(size = 10),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      # axis.title.y = element_blank(),
      # axis.text.y = element_blank()
    ) 
    # ggtitle(paste0('Patient:',p))
  # coord_flip()
  
  
  p2
  
  p3 <- ggplot()+
    # ylim(-1e18,-1e1) +
    scale_y_continuous(limits = c(-1e18,-10),expand=c(0,0))+
    scale_color_manual(values = c("Down-regulated"='blue',"Up-regulated"='red'))+
    geom_text_repel(
      data=top10_sig,
      aes(x=label,y=logFC,label=gene),
      force = 1.2,
      size=3,
      # arrow = arrow(length = unit(0.008, "npc"),
      #               type = "open", ends = "last")
    )+
    # geom_col(data = dfbar,
    #          mapping = aes(x = x,y = y),
    #          fill = "#dcdcdc",alpha = 0.6)+
    # geom_col(data = dfbar1,
    #          mapping = aes(x = x,y = y),
    #          fill = "#dcdcdc",alpha = 0.6)+
    geom_jitter(data = output_sig,
                aes(x = label, y = logFC, color = expression),
                size = 0.6,
                width = 0.3) +
    theme_minimal()+
    theme(
      axis.title = element_text(size = 13,
                                color = "black",
                                face = "bold"),
      axis.line.y.left = element_line(color = "black",
                                      size = 1.2),
      axis.line.y = element_blank(),
      # axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = c(1,0),
      legend.text = element_text(size = 10),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # axis.text.y = element_blank()
    )
  # )+
  # ylab('')+
  # coord_flip()


  p3

  p4 <- ggarrange(p1,p2,p3,heights=c(1/5, 3/5,1/5),ncol = 1, nrow = 3,common.legend = TRUE,legend="right",align = "v")
  p4
  print("p4 is or not null")
  is.null(p4)
  ggsave(paste(out_dir,paste(sample,"multiple_volcano.pdf",sep="_"),sep="/"),p4,width = 5+2*length(unique(output_sig$label)),height=15,limitsize = F,device="pdf")
  ggsave(paste(out_dir,paste(sample,"multiple_volcano.png",sep="_"),sep="/"),plot = p4,width = 5+2*length(unique(output_sig$label)),height=15,limitsize = F,device="png")
}

#p <- patient[1]
for (p in patient) {
  # Load clonetype data of the patient
  sample <- pat[which(pat$sample_title==p),"run_accession"]
  clone <- read.csv(paste(clones_dir,paste(sample,"_cells_clones_rename_barcodes_infomation.csv",sep=""),sep="/"),row.names = 1)

  # match cells to corresponding patient
  data <- subset(adata,barcodes %in% rownames(clone))
  meta_ <- subset(meta,barcodes %in% rownames(clone))
  # Add clonetype information to meta.data
  data <- AddMetaData(data,clone)
  residual_clones <- sort(unique(data$clone_id))
  output <- data.frame()
  # cl <- 0
   for (i in seq(1, length(residual_clones) - 1, 1)) {
    print(paste(residual_clones[i+1],residual_clones[i],sep="vs"))
    # select data containing two clonetype
    sub <- subset(data,clone_id %in% c(residual_clones[i],residual_clones[i+1]))
    
    # extract count matrix
    count <- as.matrix(sub@assays$RNA@data)
    
    label <- sub[['clone_id']]
    
    batch <- data.frame(barcodes=rownames(label),batch=1)
    
    df = model.matrix(~clone_id , data=label)
    
    # nebula
    # data_g = group_cell(count=count,id=batch[[2]],pred=df)
    
    re = nebula(count, batch[[2]], pred=df)
    
    # extract result
    result <- re$summary
    result <- result[,c(8,2,6)]

    colnames(result) <- c('gene','logFC','pVal')
    rownames(result) <- result$gene
    
    # BH adjusted p value
    result$adj.P <- p.adjust(result$pVal,method = 'BH')
    
    result$'-Log10(adj.P)' <- -log10(result$adj.P)
    
    result <- result[order(result$logFC,decreasing = T),]
    
    # reverse the logFC
    #result$logFC <- -(result$logFC)
    
    # define the up/down-regulated gene
    result$expression <- ifelse(result$logFC >= 1 & result$adj.P < 0.05,"Up-regulated",
                                ifelse(result$logFC <= -1 & result$adj.P < 0.05,"Down-regulated","NS."))
    
    # add a label 
    result$label <- paste(residual_clones[i+1],residual_clones[i],sep="vs")
    
    # filter signifcant differential expressed genes
    result <- result %>% 
      filter(expression != 'NS.')
    
    output <- rbind(output, result)
  } # clonetype loop
  
  # save DEGs within one patient
  print("write start")
  write.csv(output,paste(out_dir,paste(sample,'DEG_among_clonetypes.csv',sep="_"),sep="/"), quote = F)
  print("write finish")
  output_sig <- output
  output_sig$label <- factor(output_sig$label,levels = unique(output_sig$label))
  print(head(output_sig))
  if (nrow(output_sig)>1){
    print("start plot")
    plot_multiple_volcano(output_sig,out_dir)
  }
} # patient loop



