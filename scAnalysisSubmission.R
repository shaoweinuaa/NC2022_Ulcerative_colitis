setwd("~/projects/UC")
options(stringsAsFactors = F)
library(ggsci)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(yarrr)
library(tidyverse)
library(monocle3)
library(SCENIC)
library(SingleCellExperiment)
library(cowplot)
library(gridExtra)

sample.list <- list("181045B-C-0816","181045B-T-0816","181045C-C-0911","181045C-T-0911","181045D_C-1025-1", "181045D_C-1025-2",
                    "181045D_T-1025-1","181045D_T-1025-2", "181045E_C-1101-1", "181045E_C-1101-2", "181045F_C-1108-3",
                    "181045F_C-1108-4")
obj.list <- lapply(sample.list, function(x){
  data <- Read10X(data.dir = paste0("data/",x))
  data <- CreateSeuratObject(counts = data, project = x, min.cells = 3, min.features = 200)
  data[["percent.MT"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
  data <- NormalizeData(data)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  data
})
features <- SelectIntegrationFeatures(object.list = obj.list)
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
obj.combined <- IntegrateData(anchorset = obj.anchors)

DefaultAssay(obj.combined) <- "integrated"
obj.combined <- ScaleData(obj.combined, verbose = FALSE)
obj.combined <- RunPCA(obj.combined, npcs = 30, verbose = FALSE)
obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 1:30)
obj.combined <- FindNeighbors(obj.combined, reduction = "pca", dims = 1:30)
obj.combined <- FindClusters(immune.combined)

markers <- FindAllMarkers(obj.combined, only.pos = T)
markers.top20 <- markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC)

new.ident <- c("B Cell", "Epithelial","T Cell","T Cell","B Cell","T Cell","SmoothMuscle","B Cell",
               "T Cell","Macrophage","Endothelial","Mast Cell","Macrophage","Fibroblasts","SmoothMuscle",
               "NK Cell")
names(new.ident) <- levels(obj.combined)
obj.combined <- RenameIdents(obj.combined, new.ident)
pals <- piratepal(palette = "basel", trans = .2)[c(1:3,5:10)]

DefaultAssay(obj.combined) <- "RNA"
p1 <- DimPlot(obj.combined, reduction = "umap", label = T)+scale_color_manual(values = unname(pals))+NoLegend()+NoAxes()
p2 <- FeaturePlot(obj.combined, c("CD3E","CD79A","CD68","PECAM1","EPCAM","GNLY","CPA3",
                                  "DCN","SPARC"))&NoLegend()&NoAxes()

tc <- subset(obj.combined, idents = "T Cell")
tc <- ScaleData(tc)
tc <- RunPCA(tc, features = VariableFeatures(object = tc))
tc <- FindNeighbors(tc)
tc <- FindClusters(tc, resolution = 0.2)
tc <- RunUMAP(tc, dims = 1:30)
t.pals <- unname(piratepal(palette = "bugs"))
DimPlot(tc, cols = t.pals)
t.markers <- FindAllMarkers(tc, only.pos = T)
t.markers.top20 <- t.markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC)
ndx <- c("Naive CD4+ T Cell","γδ T Cell","Cytotoxic T Cell","Treg")
names(ndx) <- levels(tc)
tc <- RenameIdents(tc, ndx)
DefaultAssay(tc) <- "RNA"
tc.p1 <- DimPlot(tc, cols=t.pals, label = T)+NoLegend()+NoAxes()
tc.p2 <- FeaturePlot(tc, c("CD8A","CD4","CD3E","FOXP3","IL2RA","TRGC1","IFNG","LEF1","SELL","GZMB","CCR7","TNF"), label=F)&NoLegend()&NoAxes()

bc <- subset(obj.combined, idents = "B Cell")
bc <- ScaleData(bc)
bc <- RunPCA(bc, features = VariableFeatures(object = bc))
bc <- FindNeighbors(bc)
bc <- FindClusters(bc, resolution = 0.3)
bc <- RunUMAP(bc, dims = 1:30)
b.pals <- wesanderson::wes_palette("GrandBudapest2",n=7,type = "continuous")
DimPlot(bc, cols = b.pals, label = T)
b.markers <- FindAllMarkers(bc, only.pos = T)
b.markers.top20 <- b.markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC)
ndx <- c("CD69+ Resident Memory B1 Cell","Memory B1 Cell","Naive B Cell","GPR183+ B Cell","Breg",
         "Memory B2 Cell","GPR183+ B Cell","CD69+ Resident Memory B2 Cell","Naive B Cell",
         "CD69+ Resident Memory B1 Cell")
names(ndx) <- levels(bc)
bc <- RenameIdents(bc, ndx)
DefaultAssay(bc) <- "RNA"
bc.p1 <- DimPlot(bc, cols=b.pals, label = T)+NoLegend()+NoAxes()
bc.p2 <- FeaturePlot(bc, c("IGHA2","IGHA1","MS4A1","IGHG2","IGHM","GPR183",
                           "HSPA1A","TNF","CD69","NFKBIA","VIM","CD27"), label=F)&NoLegend()&NoAxes()

mp <- subset(obj.combined, idents = "Macrophage")
mp <- ScaleData(mp)
mp <- RunPCA(mp, features = VariableFeatures(object = mp))
mp <- FindNeighbors(mp)
mp <- FindClusters(mp, resolution = 0.2)
mp <- RunUMAP(mp, dims = 1:30)
mp.pals <- wesanderson::wes_palette("Darjeeling2",n=5,type = "discrete")[1:3]
DimPlot(mp, cols=mp.pals)
mp.markers <- FindAllMarkers(mp, only.pos = T)
mp.markers.top20 <- mp.markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC)

ndx <- c("Resident Macrophage","Infiltrating Macrophage","Resident Macrophage","Proliferating Macrophage")
names(ndx) <- levels(mp)
mp <- RenameIdents(mp, ndx)
mp.p1 <- DimPlot(mp, cols = mp.pals, label = T)+NoLegend()+NoAxes()
DefaultAssay(mp) <- "RNA"
mp.p2 <-FeaturePlot(mp, c("ITGAM","C1QA","C1QB",'CD163',"SIGLEC1","IL1B","CD68","CCR2","MKI67"), label = F)&NoAxes()&NoLegend()

sample.info <- data.frame(id=unique(obj.combined@meta.data$orig.ident), 
                          group=c("UCSC","UC","UCSC","UC","UCSC","UCSC","UC","UC","HC","HC","HC","HC"))
sample.info$group <- factor(sample.info$group, levels = c("HC","UCSC","UC"))
sample.info <- sample.info[order(sample.info$group),]
obj.combined$group <- sample.info[match(obj.combined$orig.ident, sample.info$id),"group"]
mp$group <- sample.info[match(mp$orig.ident, sample.info$id),"group"]

###### composition comparisons
library(pheatmap)
library(colorspace)

stats.list <- lapply(list(tc,bc,mp), function(x){
  stats <- data.frame(orig.ident = x$orig.ident, ident = Idents(x))
  stats <- dplyr::count(stats, orig.ident, ident )
  total <- dplyr::count(x@meta.data, orig.ident)
  stats$total <- total[match(stats$orig.ident, total$orig.ident),"n"]
  stats$freq <- stats$n/stats$total
  stats$group <- sample.info[match(stats$orig.ident, sample.info$id),"group"]
  stats
})

stats <- Reduce(rbind, stats.list)
stats <- reshape(stats[,c("orig.ident","ident","freq")],idvar = "ident", timevar = "orig.ident", direction = "wide")
rownames(stats) <- stats$ident 
stats$ident <- NULL
colnames(stats) <- gsub("^freq\\.", "", colnames(stats))
stats[is.na(stats)] <- 0
anno_col <- data.frame(group=sample.info$group, row.names = sample.info$id)
anno_row <- data.frame(CellType=c(rep("T Cell",4), rep("B Cell",7),rep("Macrophage",3)), row.names = rownames(stats))
ann_colors = list(group = setNames(wes_palette("Darjeeling1",5,"discrete")[1:3],c("UC","HC","UCSC")),
                  CellType = setNames(pals[c(1,3,5)],c("T Cell","B Cell","Macrophage")))
pheatmap(stats[,sample.info$id[sample.info$id!="181045A-1bing"]], scale='row', gaps_row = c(4,11),gaps_col = c(4,8),
         annotation_col = anno_col,annotation_row = anno_row,
         annotation_names_col=F,annotation_names_row=F,
         annotation_colors=ann_colors,cluster_rows = FALSE,
         cluster_cols = FALSE,show_colnames = FALSE,show_rownames = TRUE,border_color = "white",
         color=colorRampPalette(divergingx_hcl(7, "Geyser"))(100))

###### SOD1/2 expression analysis in Mφ
library(ggpubr)
mp.subs <- subset(mp, idents = c("Resident Macrophage", "Infiltrating Macrophage"))
mp <- AddMetaData(mp, mp@assays$RNA@data["SOD1",],"SOD1")
mp <- AddMetaData(mp, mp@assays$RNA@data["SOD2",],"SOD2")
dat <- data.frame(id=colnames(mp), ident=Idents(mp), group=mp$group, SOD1=mp$SOD1, SOD2=mp$SOD2)
dat.m <- reshape(dat, varying = 4:5, direction = "long",v.names = "gene")
dat.m$time[dat.m$time==1] <- "SOD1"
dat.m$time[dat.m$time==2] <- "SOD2"
dat.s <- dat.m %>% group_by(ident, group, time) %>% summarise(median = median(gene))

p <- ggplot(dat.m[dat.m$ident!="Proliferating Macrophage",], aes(group, gene)) + geom_boxplot(
  aes(color = ident, fill=ident), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
) +
  geom_dotplot(
    aes(fill = ident),color="white",binwidth=.2,
    binaxis='y', stackdir='center', dotsize = 0.8,
    position = position_dodge(0.8)
  )+
  geom_line(data=dat.s[dat.s$ident!="Proliferating Macrophage",], 
            aes(x=group, y=median,group = ident, color=ident), linetype=4, size=2)+
  scale_fill_manual(values = mp.pals)+
  scale_color_manual(values = mp.pals)+
  stat_compare_means(aes(group = ident), label = "p.signif")+
  facet_wrap(~time,ncol = 1)+theme_classic()+theme(legend.title = element_blank())+xlab("")+
  ylab("Expression")

###### pathway activities
library(GSVA)
library(ggridges)
library(usethis) 
usethis::edit_r_environ()
genesets <- list()
genesets[["ROS"]] <- read.table("geneset-2.txt", skip = 2)$V1
test <- apply(obj.combined@assays$RNA@data, 1, function(x)sum(x>0))
test <- names(test)[test>=5000]
obj.list <- SplitObject(obj.combined, split.by = "ident")
gsea.lists <- lapply(obj.list, function(x){
  gsva(as.matrix(x@assays$RNA@data[test,]), genesets, min.sz=1, max.sz=500)
})

gsea <- c()
for(x in names(gsea.lists)){
  gsea <- c(gsea,setNames(gsea.lists[[x]], colnames(gsea.lists[[x]])))
}

obj.combined <- AddMetaData(obj.combined, gsea, "ROS")
obj.combined$group.subset <- as.character(obj.combined$Idents)
obj.combined$group.subset[names(Idents(tc))] <-as.character(Idents(tc))
obj.combined$group.subset[names(Idents(bc))] <-as.character(Idents(bc))
obj.combined$group.subset[names(Idents(mp))] <-as.character(Idents(mp))
obj.combined$group.detail <- paste0(obj.combined$group.subset,"(",obj.combined$group,")")

pals <- wes_palette("Darjeeling1",5,"discrete")[c(c(2,3,1))]
for(item in unique(obj.combined$group.subset)){
  dat <- obj.combined@meta.data[obj.combined$group.subset== item,]
  dat$group.detail <- factor(dat$group.detail, 
                             levels = paste0(item,"(",c("UC","UCSC","HC"),")"))
  p <- ggplot(dat, aes(x = ROS, y = group.detail, fill=group, alpha=.9)) + 
    geom_density_ridges(scale = 2, quantile_lines=TRUE,color="white",
                        quantile_fun=function(x,...)median(x))+ylab("")+xlab("GSEA Score(ROS)")+
    theme_classic()+scale_fill_manual(values = pals)+theme(legend.position = "none")
}

###### Differential expression analysis
library(ggrepel)
mp.uc <- subset(mp.subs, group=="UC")
mp.hc <- subset(mp.subs, group=="HC")
mp.ucsc <- subset(mp.subs, group=="UCSC")

mp.if <- subset(mp.subs, idents = "Infiltrating Macrophage")
mp.rs <- subset(mp.subs, idents = "Resident Macrophage")

Idents(mp.if) <- "group"
Idents(mp.rs) <- "group"

diffExpr <- FindMarkers(mp.uc, ident.1 = "Infiltrating Macrophage", ident.2 = "Resident Macrophage",
                        logfc.threshold=0, min.pct = .5)

diffExpr$direction <- "NotSig"
diffExpr[diffExpr$p_val_adj<0.05 & diffExpr$avg_log2FC< -log2(1.5),"direction"] <- "Down"
diffExpr[diffExpr$p_val_adj<0.05 & diffExpr$avg_log2FC> log2(1.5),"direction"] <- "Up"
table(diffExpr$direction)
diffExpr$symbol <- NA
diffExpr$symbol[abs(diffExpr$avg_log2FC)>log2(1.5) & diffExpr$p_val_adj<0.05] <- 
  rownames(diffExpr[abs(diffExpr$avg_log2FC)>log2(1.5) & diffExpr$p_val_adj<0.05,])

pals <- wes_palette("Zissou1", 2, type = "continuous")
p <- ggplot(data=diffExpr, aes(x=avg_log2FC, y=-log10(p_val_adj), col=direction, label=symbol)) + 
  geom_point() + 
  theme_classic() +
  geom_text_repel(max.overlaps=20) +
  scale_color_manual(values=c(pals[1], "grey80", pals[2])) +
  geom_vline(xintercept=c(-log2(1.5), log2(1.5)), col="grey50",linetype = 'dotted') +
  geom_hline(yintercept=-log10(0.05), col="grey50",linetype = 'dotted')+
  theme(legend.position = "none")

library(clusterProfiler)
library(msigdbr)
h_df <- msigdbr(species = "Homo sapiens")
head(h_df, 2) %>% as.data.frame
h.hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

symbols.up <- rownames(diffExpr[diffExpr$avg_log2FC>log2(1.5) & diffExpr$p_val_adj<0.05,])
symbols.down <- rownames(diffExpr[diffExpr$avg_log2FC< -log2(1.5) & diffExpr$p_val_adj<0.05,])

h.up <- enricher(symbols.up, TERM2GENE=h.hallmark)
h.down <- enricher(symbols.down, TERM2GENE=h.hallmark)

h.up <- data.frame(group="UpRegulated", h.up@result[h.up@result$p.adjust<.05, c("ID","p.adjust","geneID")])
h.down <- data.frame(group="DownRegulated", h.down@result[h.down@result$p.adjust<.05, c("ID","p.adjust","geneID")])

diff.top <- rbind(h.up, rev(h.down))
diff.top$ID <- gsub("^HALLMARK_","",diff.top$ID )
diff.top$ID <- gsub("_"," ",diff.top$ID )

diff.top$ID <- factor(diff.top$ID, levels = rev(diff.top$ID))
diff.top$x <- -log10(diff.top$p.adjust)
diff.top$count <- str_count(diff.top$geneID,"/")+1

library(stringi)
library(stringr)
p <- ggplot(diff.top,aes(x=ID,y=x))+geom_bar(stat = "identity", aes(fill=group))+coord_flip()+
  xlab("")+ylab("-log10 (adjusted p-value)")+theme_classic()+ylim(0,37)+
  scale_fill_manual(values = pals,name="")+  
  geom_text(aes(label = count, x = ID, y = x),hjust = 1.2,color="white")

######  SCENIC analysis
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
}

library(SCENIC)
exprMat <- as.matrix(mp@assays$RNA@counts)
cellInfo <- data.frame(seuratCluster=Idents(mp), group=paste0(Idents(mp),"-",mp$group))
scenicOptions <- initializeScenic(org="hgnc", dbDir=".", nCores=10)
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

exprMat_log <- log2(exprMat+1)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings
saveRDS(scenicOptions, file="scenicOptions.Rds") 

motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox8"]
viewMotifs(tableSubset) 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "group"])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
regulonTargetsInfo[regulonTargetsInfo$gene=="TNF",]
regulonTargetsInfo[regulonTargetsInfo$gene=="IL1B",]

nCores <- getSettings(scenicOptions, "nCores")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
thresholds <- loadInt(scenicOptions, "aucell_thresholds")
thresholds <- getThresholdSelected(thresholds)
regulonsCells <- setNames(lapply(names(thresholds), 
                                 function(x) {
                                   trh <- thresholds[x]
                                   names(which(getAUC(regulonAUC)[x,]>trh))
                                 }),names(thresholds))

regulonActivity <- reshape2::melt(regulonsCells)
binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
class(binaryRegulonActivity) <- "matrix"
binaryRegulonActivity_nonDupl <- binaryRegulonActivity[which(rownames(binaryRegulonActivity) %in% onlyNonDuplicatedExtended(rownames(binaryRegulonActivity))),]
AUC_values <- getAUC(regulonAUC)
auc_assay <- CreateAssayObject(counts = AUC_values[,rownames(mp.subs@meta.data)])
mp.subs[["AUC"]] <- auc_assay
Assays(mp.subs)
DefaultAssay(mp.subs) <- "AUC"
mp.subs <- AddMetaData(mp.subs, cellInfo[rownames(mp.subs@meta.data),"group"], "group")
Idents(mp.subs) <- "group"
auc.markers <- FindAllMarkers(mp.subs, only.pos = T,logfc.threshold = .01, min.pct = 0.05)
auc.markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC) -> best.tfs
AUC_values <-AUC_values[unique(gsub("-","_",best.tfs$gene)),
                        rownames(mp.subs@meta.data)[order(Idents(mp.subs))]]

anno_col <- cellInfo[,"group", drop=F]
ann_colors = list(group = setNames(brewer.pal(6,"Set3"),
                                   c("Infiltrating Macrophage-HC","Infiltrating Macrophage-UCSC",
                                     "Infiltrating Macrophage-UC","Resident Macrophage-HC",
                                     "Resident Macrophage-UCSC","Resident Macrophage-UC")))
pheatmap(AUC_values, annotation_col = anno_col, show_colnames = FALSE, 
         annotation_names_col = F, scale = "row", cluster_rows = T, 
         cluster_cols = T, annotation_colors = ann_colors)

##### using Nitchnet for cell-cell interaction analysis
ccg.list <- read.table("GeneList.txt", header = T, sep="\t")
table(ccg.list$Category)
ccg.list <- ccg.list[grep("Chemokine|Cytokine",ccg.list$Category), ]

library(nichenetr)
ligand_target_matrix = readRDS("~/Downloads/ligand_target_matrix.rds")
lr_network = readRDS("~/Downloads/lr_network.rds")
weighted_networks = readRDS("~/Downloads/weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

obj.combined$Idents <- Idents(obj.combined)
Idents(obj.combined) <- "group.subset"
obj.uc <- subset(obj.combined, group=="UC")
obj.hc <- subset(obj.combined, group=="HC")
obj.ucsc <- subset(obj.combined, group=="UCSC")

##### Analysis for Tcells with Mφ
receiver = levels(obj.uc)[c(1,5,11,19)]
expressed_genes_receiver_uc = receiver %>% unique() %>% lapply(get_expressed_genes, obj.uc, 0.03, assay_oi="RNA")
expressed_genes_receiver_uc = expressed_genes_receiver_uc %>% unlist() %>% unique()
background_expressed_genes = expressed_genes_receiver_uc %>% .[. %in% rownames(ligand_target_matrix)]

sender = levels(obj.uc)[c(14,18)]
expressed_genes_sender_uc = sender %>% unique() %>% lapply(get_expressed_genes, obj.uc, 0.5, assay_oi="RNA") 
expressed_genes_sender_uc = expressed_genes_sender_uc %>% unlist() %>% unique()

geneset_oi = Reduce(intersect, list(expressed_genes_receiver_uc, ccg.list$Symbol, rownames(ligand_target_matrix)))
ligands = lr_network %>% pull(from) %>% unique() %>% intersect(ccg.list$Symbol)
receptors = lr_network %>% pull(to) %>% unique() %>% intersect(ccg.list$Symbol)
expressed_ligands = intersect(ligands,expressed_genes_sender_uc)
expressed_receptors = intersect(receptors,expressed_genes_receiver_uc)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,
                                                                 geneset = geneset_oi, 
                                                                 ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows() %>% drop_na()


##### Analysis for Bcells with Mφ
receiver = levels(obj.uc)[c(4,7,13,15,17,20,22)]
expressed_genes_receiver_uc = receiver %>% unique() %>% lapply(get_expressed_genes, obj.uc, 0.03, assay_oi="RNA")
expressed_genes_receiver_uc = expressed_genes_receiver_uc %>% unlist() %>% unique()
background_expressed_genes = expressed_genes_receiver_uc %>% .[. %in% rownames(ligand_target_matrix)]

sender = levels(obj.uc)[c(14,18)]
expressed_genes_sender_uc = sender %>% unique() %>% lapply(get_expressed_genes, obj.uc, 0.5, assay_oi="RNA") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender_uc = expressed_genes_sender_uc %>% unlist() %>% unique()

geneset_oi = Reduce(intersect, list(expressed_genes_receiver_uc, ccg.list$Symbol, rownames(ligand_target_matrix)))
ligands = lr_network %>% pull(from) %>% unique() %>% intersect(ccg.list$Symbol)
receptors = lr_network %>% pull(to) %>% unique() %>% intersect(ccg.list$Symbol)
expressed_ligands = intersect(ligands,expressed_genes_sender_uc)
expressed_receptors = intersect(receptors,expressed_genes_receiver_uc)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,
                                                                 geneset = geneset_oi, 
                                                                 ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows() %>% drop_na()
