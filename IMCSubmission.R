setwd("~/projects/UC/")
options(stringsAsFactors = F)
library(Rphenograph)
library(doParallel)
library(reshape)
library(ggplot2)
library(Rtsne)
library(pheatmap)
library(ggsci)
library(wesanderson)
library(dplyr)
library(gridExtra)
library(ggvoronoi)
library(harmony)
library(imcRtools)
library(batchelor)
library(batchelor)
library(ggrepel)


files <- list.files("NC/NC/csvFile/", "*csv", full.names = T)
xys = foreach(file=files, .combine=rbind) %dopar%
  {
    df <- read.table(file, header = T, sep=",")
    df[,1:43]
  }
rownames(xys) <- paste(xys$ROI, xys$CellId, sep="-")
xys[,5:43] <- apply(xys[,5:43], 2, function(x){
  x <- (x-quantile(x,.01))/(quantile(x,.99)-quantile(x,.01))
})

set.seed(42)
harmony_emb <- HarmonyMatrix(xys[,5:43], ifelse(grepl("gut", rownames(xys)),"A","B"),
                             do_pca = TRUE)
Rphenograph_out <- Rphenograph(harmony_emb, k=100)
xys$Phenograph <- factor(paste0("Cluster",membership(Rphenograph_out[[2]])))
tsne_out <- Rtsne(harmony_emb,check_duplicates=FALSE, set.seed=42)

anno.tbl <- data.frame(id = paste0("Cluster", c(seq(1,22))), 
                       annotation = c("ResidentMacrophages","Vimentin+MesenchymalCells","ResidentMacrophages","Gata3+TCells",
                                      "InfiltratingMacrophages/TNFα+TCells","ExhausticTCells","Neutrophils",
                                      "Vimentin+MesenchymalCells","InfiltratingMacrophages","Vimentin+MesenchymalCells",
                                      "Epithelials/ImmuneSuppressiveCells","Gata3+CD8+TCells",
                                      "InfiltratingMacrophages/TNFα+TCells/NK","Vimentin+MesenchymalCells","Lineage-Cells",
                                      "Tregs","NK/BCells","Vimentin+MesenchymalCells","TNFα+Epithelials",
                                      "Fibroblasts/Endothelials",
                                      "ResidentMacrophages/Epithelials","Epithelials"))
name.tbl <- data.frame(ref=colnames(xys)[5:43],alt=c(gsub(".*_(.*)","\\1", colnames(xys)[5:12]),
                                                     "Caspase3","CD31","E-Cadherin","B7-H4","VISTA",
                                                     gsub(".*_(.*)","\\1", colnames(xys)[18:20]),"PD-L1",
                                                     "LAG3","CD68","CD11b","CD20","CD11c","CD15","Granzyme B","PD1",
                                                     "Ki67","GATA3","HLA-DR","CD45RA","CD3","TNFα","IL-1β","CD45RO","CD57",
                                                     "CD25","Pan-Keratin","αSMA","Vimentin","CD45"))
xys$Phenograph <- anno.tbl[match(xys$Phenograph,anno.tbl$id),"annotation"]
pals <- setNames(pal_d3("category20")(17),sort(unique(xys$Phenograph)))
exprs.m <- aggregate(xys[,5:43], by=list(xys$Phenograph), mean)
rownames(exprs.m) <- exprs.m$Group.1
exprs.m$Group.1 <- NULL
exprs.m$Phenograph <- NULL
colnames(exprs.m) <- name.tbl[match(colnames(exprs.m), name.tbl$ref),"alt"]
hm <- pheatmap(exprs.m, border_color = "white", scale = "none",
               color= wes_palette("Zissou1", 100, type = "continuous"))

dat <- data.frame(tsne_out$Y)
dat$Phenograph <- factor(xys$Phenograph)
colnames(dat)[1:2] <- c("TSNE1","TSNE2")
pos.df <- aggregate(dat[,c("TSNE1","TSNE2")], median, by=list(dat$Phenograph))
pos.df$key <- paste(seq(1,17), pos.df$Group.1, sep=":")
pos.df$label <- seq(1,17)
pals <- setNames(pal_d3("category20")(17),sort(unique(xys$Phenograph)))
p <- ggplot(data = dat, aes(TSNE1,TSNE2))+geom_point(aes(color=Phenograph))+ 
  geom_text(data = pos.df, aes(x=TSNE1,y=TSNE2,label=label),size=8,color="grey10")+
  scale_color_manual(values = pals, name="", labels=pos.df$key)+
  theme_classic()+theme(legend.text = element_text(size=14))
cbind(xys[,5:43],dat) %>% select(-Phenograph) %>% melt(.,id.vars=c("TSNE1","TSNE2")) -> dat.draw
dat.draw$variable <- name.tbl[match(dat.draw$variable, name.tbl$ref),"alt"]
p.list <- lapply(unique(dat.draw$variable), function(i) {
  ggplot(data = dat.draw[dat.draw$variable==i,], aes(TSNE1,TSNE2))+
    geom_point(aes(color=value),size=.6)+ 
    scale_color_gradientn(colours = hcl.colors(12, "Geyser"))+ggtitle(i)+
    theme_void()+ theme(plot.background = element_rect(fill = "white",colour = NA),
                        panel.background = element_rect(fill = "white",colour = NA),
                        plot.title = element_text(hjust = 0.5, face="bold", size=15),
                        legend.position = "none")
})
p <- do.call(grid.arrange, c(p.list, ncol=10))

hm.orders <- rownames(exprs.m)[hm$tree_row$order]
ct.stats <- data.frame(roi=gsub("(.*)-.*","\\1",rownames(xys)), Phenograph=xys$Phenograph)
roi.total <- dplyr::count(ct.stats, roi)
ct.stats <- dplyr::count(ct.stats, roi, Phenograph)
ct.stats$total <- roi.total[match(ct.stats$roi, roi.total$roi),"n"]
ct.stats$freq <- ct.stats$n/ct.stats$total
ct.stats <- reshape(ct.stats[,c(1,2,5)],idvar = "roi", timevar = "Phenograph", direction = "wide")
ct.stats[is.na(ct.stats)] <- 0
ct.stats <- melt(ct.stats)
ct.stats$variable <- gsub("freq\\.","", ct.stats$variable)
ct.stats$variable <- factor(ct.stats$variable, levels = hm.orders)
ggplot(ct.stats, aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill=variable),color="white") + 
  geom_jitter(aes(fill=variable),color="white", pch=21) + theme_classic()+
  scale_fill_manual(values = pals[hm.orders])+
  theme(legend.position = "none",axis.text.x = element_blank())+
  ylab("Frequency")+xlab("")


library(factoextra)
stats <- data.frame(ct.stats[,c("roi","variable", "value")])
stats.w <- reshape(stats, timevar = "variable", idvar = "roi", direction = "wide")
stats.w[is.na(stats.w)] <- 0
df <- scale(stats.w[,2:18],center = TRUE,scale = TRUE)
df <- data.frame(df)
rownames(df) <- stats.w$roi
colnames(df) <- as.character(unique(ct.stats$Phenograph))

res.dist = dist(x = df,method = "euclidean")
res.hc <- hclust(d = res.dist,method = "ward.D")

library(dendextend)
my_colors <- ifelse(grepl("normal|_3_", rownames(df)),wes_palette("Darjeeling1",5,"discrete")[2],
                    wes_palette("Darjeeling1",5,"discrete")[1])
as.dendrogram(res.hc) %>%
  set("labels_cex", 0)  %>%
  set("labels", "")  %>%
  set("leaves_pch", 19)  %>%
  set("nodes_cex", 0.7) %>% 
  set("branches_lwd",2) %>%
  plot(axes=FALSE)
colored_bars(colors = my_colors, dend = as.dendrogram(res.hc))
library("FactoMineR")
res.pca <- PCA(df,  graph = FALSE)
p  <- fviz_pca_ind(res.pca,label = "none", 
                   geom="point",pointsize = 2,
                   fill.ind = as.factor(ifelse(grepl("normal|_3_", rownames(df)),"Normal","UC")), 
                   col.ind = "white",
                   pointshape = 21, 
                   palette = c(wes_palette("Darjeeling1",5,"discrete")[2],wes_palette("Darjeeling1",5,"discrete")[1]),
)

library(imcRtools)
nbrs <- data.frame(from=as.character(), to=as.character())
sa <- data.frame(group_by=as.character(),from_label=as.character(),to_label=as.character(),ct=as.numeric(),
                 p_gt=as.numeric(),p_lt=as.numeric(),interaction=as.logical(),p=as.numeric(),sig=as.logical(),
                 sigval=as.numeric())   
for (file in files){  
  dat <- read.table(file, header = T, sep=",")
  rownames(dat) <- paste(dat$ROI, dat$CellId, sep="-")
  cd <- DataFrame(x = dat$X_position, y = dat$Y_position, z=rownames(dat), 
                  ImageNb=gsub("NC\\/NC\\/csvFile\\/\\/(.*)\\.csv","\\1", file),
                  Pos_X=dat$X_position, Pos_Y= dat$Y_position, 
                  CellType=xys[match(rownames(dat), rownames(xys)),"Phenograph"])
  spe1 <- SpatialExperiment(
    assay = t(dat[,5:43]), 
    colData = cd, 
    spatialDataNames = "z", 
    spatialCoordsNames = c("Pos_X", "Pos_Y"))
  
  spe1 <- buildSpatialGraph(spe1, img_id = "ImageNb",type = "knn",k = 20)
  tmp <- data.frame(colPair(spe1, "knn_interaction_graph"))
  tmp <- apply(tmp, 2, function(x)
    paste(gsub("NC\\/NC\\/csvFile\\/\\/(.*)\\.csv","\\1", file),x, sep="-"))
  nbrs <- rbind(nbrs, tmp)
  out <- testInteractions(spe1,
                          group_by = "ImageNb",
                          label = "CellType",
                          method = "classic",
                          colPairName = "knn_interaction_graph")
  sa <- rbind(sa, data.frame(out))
}

nbrs.df <- dplyr::count(nbrs, from, Phenograph)
nbrs.df <- reshape(nbrs.df, idvar = "from", timevar = "Phenograph", direction = "wide")
rownames(nbrs.df) <- nbrs.df$from
nbrs.df$from <- NULL
nbrs.df[is.na(nbrs.df)] <-0
set.seed(42)
km.res <- kmeans(scale(nbrs.df), 15, nstart = 42)
nbrs.df$Phenograph <-paste0("CN",km.res$cluster)

nbrs.m <- aggregate(nbrs.df, by=list(nbrs.df$Phenograph), mean)
rownames(nbrs.m) <- nbrs.m$Group.1
nbrs.m$Group.1 <- NULL
nbrs.m$Phenograph <- NULL
colnames(nbrs.m) <- gsub("n\\.","",colnames(nbrs.m))
hm <- pheatmap(nbrs.m, border_color = "white", scale = "column",
               color= wes_palette("Zissou1", 100, type = "continuous"))


library(ggvoronoi)
pals <- setNames(pal_d3("category20")(15),paste0("CN",1:15))
vrns <- cbind(nbrs.df$Phenograph, xys[match(rownames(nbrs.df), rownames(xys)),1:4])
colnames(vrns)[1] <- "Phenograph"
rois <- unique( vrns$ROI)
vrns$Phenograph <- factor(vrns$Phenograph, levels = paste0("CN",1:15))
for (roi in rois){
  ggplot(vrns[vrns$ROI==roi,]) +
    geom_point(aes(X_position,Y_position, color=Phenograph), alpha=1)+
    geom_voronoi(alpha=.8, aes(X_position,Y_position,fill=Phenograph))+theme_classic()+
    scale_fill_manual(values = pals,name="")+scale_color_manual(values = pals, guide="none")+
    coord_flip()+
    scale_x_reverse()+theme_void()
}

##### spatial interaction analysis
cor.mat <- stats.w[,2:18]
cor.order <- rev(unique(xys$Phenograph)[c(1,16,17,8,14,3,5,9,12,13,4,7,10,6,2,15,11)])

rownames(cor.mat) <- stats.w$roi
colnames(cor.mat) <- gsub("value\\.","", colnames(cor.mat))

corn.mat <- cor.mat[grep("normal|_3_", rownames(cor.mat),invert = T),]
corn.mat <- cor(corn.mat, method = "pearson", use = "complete.obs")
corn.mat[is.na(corn.mat)] <- 0

corn.mat <- melt(corn.mat)
corn.mat$X1 <- factor(corn.mat$X1, levels = cor.order)
corn.mat$X2 <- factor(corn.mat$X2, levels = cor.order)
colnames(corn.mat)[1:2] <- c("from_label","to_label")

san <- sa[grep("normal|_3_", sa$group_by, invert = T),]
nt <- length(unique(san$group_by))
san <- san[san$sig==TRUE,]
attr <- dplyr::count(san[san$sigval==1,],from_label,to_label)
avod <- dplyr::count(san[san$sigval==-1,],from_label,to_label)
corn.mat <- left_join(corn.mat, attr, by=c("from_label", "to_label"))
corn.mat <- left_join(corn.mat, avod, by=c("from_label", "to_label"))
corn.mat[is.na(corn.mat)] <- 0
corn.mat$group <- ifelse(corn.mat$n.x>corn.mat$n.y,"Interaction","Avoidance")
corn.mat[corn.mat$n.x==corn.mat$n.y,"group"] <- 
  ifelse(corn.mat[corn.mat$n.x==corn.mat$n.y,"value"] >0,"Interaction","Avoidance")
corn.mat$percent <- apply(corn.mat[,4:5],1,function(x)max(x)/nt) 

p <- ggplot(corn.mat, aes(x = from_label, y = to_label, fill = value)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  geom_point(aes(x = from_label, y = to_label, color = group, size=ifelse(percent==0,NA,percent)))+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="bold",size=10),
        axis.text.y = element_text(face="bold",size=10),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        panel.border = element_blank())+
  scale_fill_gradientn(colors=hcl.colors(12, "Fall"),limits=c(-1,1), name="Correlation Coefficients")+
  scale_color_manual(values=rev(wes_palette("Darjeeling1")), name="")+
  scale_size(name="%ROI with significant\nInteraction/Avoidance",limits = c(0,1))+
  xlab("")+ylab("")+scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))



