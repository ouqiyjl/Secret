#################结果二 伪时间轨迹分析
rm(list = ls())
#################加载R包
library(monocle3)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(ggplot2)
library(tidyverse)
library(data.table)
library(topGO)
library(enrichplot)
#devtools::install_github('ievaKer/aPEAR')
library(aPEAR)
#################首读入数据------
cancer_name <- "ESCC"

dir.create(paste0("result2/",cancer_name))

load(paste0("result1/",cancer_name,'/',"PCD_ANN.RData"))
malignant_gene <- read_xlsx("gene_publication.xlsx")
malignant_gene <- unique(malignant_gene$gene)
malignant_gene <- malignant_gene[-grep("hsa",x = malignant_gene)]

new_data <-  AddModuleScore(object =new_data,features =list(malignant_gene) ,
                            name = "malignant_Score")
new_data <- AddMetaData(new_data,metadata = Idents(new_data),col.name = "Death_stats")

mean_score <- aggregate(new_data$malignant_Score1, by=list(type=new_data$Death_stats),mean)

###############对其进行CDS对象构建
data <- GetAssayData(new_data, assay = 'RNA', slot = 'counts')
cell_metadata <-new_data@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cell_metadata$cell_type = new_data@active.ident
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
gc()
#############monocle3 的标准化流程
cds <- preprocess_cds(cds, num_dim = 60)
############CDS 的降维
cds <- reduce_dimension(cds,preprocess_method = "PCA", reduction_method="UMAP") #preprocess_method默认是PCA
############从seurat导入整合过的umap坐标

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(new_data, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
#############进行轨迹分析
#进行聚类

cds <- cluster_cells(cds,reduction_method = c("UMAP"))
##进行为时间分析
cds <- learn_graph(cds,verbose = T,use_partition = F)

earliest_state <- mean_score$type[which.min(mean_score$x)]


get_earliest_principal_node <- function(cds, cell_type){
  cell_ids <- which(colData(cds)[, "cell_type"] == cell_type)
  
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

cds<-order_cells(cds,root_pr_nodes=get_earliest_principal_node(cds,earliest_state))
save(cds,file = paste0("result2/",cancer_name,"/CDS.RData"))
###########画伪时间轨迹
pseudotime_plot_file  <-  paste("result2",cancer_name,"pseudotime_plot.pdf",sep = "/")
pdf(file = pseudotime_plot_file,width = 11.69,height = 8.27)
plot_cells(cds, color_cells_by = "pseudotime", 
           label_cell_groups = F,
           label_leaves = T,  label_branch_points = T,
           group_label_size = 4,
           labels_per_group = 4, graph_label_size = 2,
           label_groups_by_cluster =F,
           label_principal_points =F,
           cell_size = 1.5,trajectory_graph_segment_size = 1,
           reduction_method = "UMAP",trajectory_graph_color="black")+
  scale_color_gradient("Pseudotime",
                       low = "#F0F921FF",
                       high = "#8707A6FF")+
  theme(axis.title.x=element_text(size=15,color = "black"),
        axis.title.y=element_text(size=15,color = "black"),
        axis.text.y = element_text(size=16,color = "black"),
        axis.text.x = element_text(size=16,color = "black"))+
  theme(legend.text = element_text(size =12,color = "black"),
        legend.title =element_text(size =12,color = "black") )
dev.off()
##########画恶化的打分加轨迹
malignant_plot_file  <-  paste("result2",cancer_name,"malignant_plot.pdf",sep = "/")
pdf(file = malignant_plot_file,width = 11.69,height = 8.27)
plot_cells(cds, color_cells_by = "malignant_Score1", 
           label_cell_groups = F,
           label_leaves = T,  label_branch_points = T,
           group_label_size = 4,
           labels_per_group = 4, graph_label_size = 2,
           label_groups_by_cluster =F,
           label_principal_points =F,
           cell_size = 1.5,trajectory_graph_segment_size = 1,
           reduction_method = "UMAP",trajectory_graph_color="#482936")+
  scale_color_gradient("Malignant_Score",
                       low = "#FFFFFF",
                       high ="#00008B")+
  theme(axis.title.x=element_text(size=15,color = "black"),
        axis.title.y=element_text(size=15,color = "black"),
        axis.text.y = element_text(size=16,color = "black"),
        axis.text.x = element_text(size=16,color = "black"))+
  theme(legend.text = element_text(size =12,color = "black"),
        legend.title =element_text(size =12,color = "black"))
dev.off()


######死亡状态
color_PCD <- c("Alkaliptosis"="#41AB5D","Anoikis" = "#A35D34","Autophagy"="#ABCADE",
               "Cuproptosis"="#8DD3C7",
               "Entotie_cell_death"="#B3DE69","Extrinsic_apoptosis"="#3B77AC",
               "Ferroptosis"="#BEBADA",
               "Immunogenic_cell_death"="#C6DBEF","Intrinsic_apoptosis"="#80B1D3",
               'Lysosome_dependent_cell_death'="#c06f98",
               "Necroptosis"="#FB8072","Netotic_cell_death"="#EF9BB5",
               "Oxeiptosis"="#F3B161",
               "Parthanatos"="#B7DA90","Pyroptosis" = "#EB9B98",
               "Extrinsic_Intrinsic_apoptosis"="#CF352B",
               "None"="#8A959B")

pcd_plot_file  <-  paste("result2",cancer_name,"pcd_plot.pdf",sep = "/")
colnames(colData(cds))[which(colnames(colData(cds))=="cell_type")] <- "Death_state"
pdf(file = pcd_plot_file,width = 11.69,height = 8.27)
plot_cells(cds, color_cells_by = "Death_state", 
           label_cell_groups = F,
           label_leaves = T,  label_branch_points = T,
           group_label_size = 4,
           labels_per_group = 4, graph_label_size = 2,
           label_groups_by_cluster =F,
           label_principal_points =F,
           cell_size = 1.5,trajectory_graph_segment_size = 1,
           reduction_method = "UMAP",trajectory_graph_color="#482936")+
  theme(axis.title.x=element_text(size=15,color = "black"),
        axis.title.y=element_text(size=15,color = "black"),
        axis.text.y = element_text(size=16,color = "black"),
        axis.text.x = element_text(size=16,color = "black"))+
  theme(legend.text = element_text(size =12,color = "black"),
        legend.title =element_text(size =12,color = "black"))+
  scale_color_manual(values = color_PCD)
dev.off()
##################计算轨迹差异基因

pr_graph_test_res <- graph_test(cds, neighbor_graph="principal_graph",cores = 8)
save(pr_graph_test_res,file = paste0("result2/",cancer_name,"/pr_graph_test_res.RData"))
##################提取差异基因
pr_graph_test_res_1 <- pr_graph_test_res[- which(abs(pr_graph_test_res$morans_I)<0.1),]
pr_deg_ids <- row.names(subset(pr_graph_test_res_1, q_value < 0.001))
#########计算亚群之间的基因共表达模块
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)),verbose = T)
save(gene_module_df,file =paste0("result2/",cancer_name,"/gene_module.RData") )
#########绘制gene module
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$Death_state)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

annotation_col = data.frame(
  cell_death  =names(color_PCD))
rownames(annotation_col) = names(color_PCD)

ann_colors = list(
  cell_death =  c("Alkaliptosis"="#41AB5D","Anoikis" = "#A35D34","Autophagy"="#ABCADE",
                  "Cuproptosis"="#8DD3C7",
                  "Entotie_cell_death"="#B3DE69","Extrinsic_apoptosis"="#3B77AC",
                  "Ferroptosis"="#BEBADA",
                  "Immunogenic_cell_death"="#C6DBEF","Intrinsic_apoptosis"="#80B1D3",
                  'Lysosome_dependent_cell_death'="#c06f98",
                  "Necroptosis"="#FB8072","Netotic_cell_death"="#EF9BB5",
                  "Oxeiptosis"="#F3B161",
                  "Parthanatos"="#B7DA90","Pyroptosis" = "#EB9B98",
                  "Extrinsic_Intrinsic_apoptosis"="#CF352B",
                  "None"="#8A959B"))


p1 <- pheatmap::pheatmap(agg_mat,
                         scale="column", 
                         clustering_method="ward.D2",
                         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = annotation_col,annotation_colors = ann_colors)

dev.new()
pcd_plot_file  <-  paste("result2",cancer_name,"modules_plot.pdf",sep = "/")
pdf(file = pcd_plot_file,width = 8.27,height = 11.69)
print(p1)
dev.off()
##################绘制基因的GO——BP
enrich <- enrichGO(pr_deg_ids, OrgDb = org.Hs.eg.db,
                   ont = 'BP',keyType = "SYMBOL",
                   pvalueCutoff  = 0.01)

save(enrich,file = paste("result2",cancer_name,"GO_BP.RData",sep = "/"))
########### 画网络图
enrichmentData <- enrich@result %>%
  as.data.table() %>%
  .[ , list(Description, pathwayGenes = geneID, p.adjust, Size = Count) ]

enrichmentData <- enrichmentData[order(enrichmentData$p.adjust)[1:300],]

enrichmentData <- as.data.frame(enrichmentData)
enrichmentData$p.adjust <- -log10(enrichmentData$p.adjust )
p1 <- enrichmentNetwork(enrichmentData, colorBy = 'p.adjust',
                        colorType = 'nes', nodeSize = 'Size', 
                  fontSize = 3.5,
                  verbose = TRUE,
                  repelLabels = TRUE,
                  drawEllipses = TRUE)+
  scale_color_gradient2("-log10(p.adjust)",
                        high = "#B83D3D",
                        mid="white",
                        low='#1A5592')

pcd_plot_file  <-  paste("result2",cancer_name,"GO_BP_plot_1.pdf",sep = "/")
pdf(file = pcd_plot_file,width = 11.69,height = 8.27)
print(p1)
dev.off()


#########网络图2
edo <- pairwise_termsim(enrich)
p2 <- emapplot(edo,layout="kk",
         showCategory = 30)

pcd_plot_file  <-  paste("result2",cancer_name,"GO_BP_plot_2.pdf",sep = "/")
pdf(file = pcd_plot_file,width = 11.69,height = 8.27)
print(p2)
dev.off()

###########画BP 柱状图
enrichmentData <- enrich@result %>%
  as.data.table() %>%
  .[ , list(Description, pathwayGenes = geneID, p.adjust,  Count) ]

enrichmentData <- enrichmentData[order(enrichmentData$p.adjust)[1:300],]

enrichmentData <- as.data.frame(enrichmentData)
############
BP <- enrichmentData[1:20,]


#指定绘图顺序（转换为因子）：
BP$term <- factor(BP$Description, levels = rev(BP$Description))


#GO富集柱形图绘制

#自定义主题：
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))

#创建一个可套用的绘图函数减少重复代码：
GO_bar <- function(x){
  y <- get(x)
  ggplot(data = y, aes(x = Count, y = term, fill = -log10(p.adjust))) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_y_discrete(labels = function(y) 
      str_wrap(y, width = 30) ) + #label换行，部分term描述太长
    labs(x = "Gene Number", y = "",
         title = paste0(x, " of GO enrichment barplot")) +
    theme_bw() +
    mytheme+
    theme(axis.title.x=element_text(color = "black"),
          axis.title.y=element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.text.x = element_text(color = "black"))+
    theme(legend.text = element_text(color = "black"),
          legend.title =element_text(color = "black") )
}

p3 <- GO_bar("BP")+scale_fill_distiller(palette = "Oranges",direction = 1)

########保存
pcd_plot_file  <-  paste("result2",cancer_name,"GO_BP_plot_3.pdf",sep = "/")
pdf(file = pcd_plot_file,width = 8.27,height =8.27 )
print(p3)
dev.off()



save.image(file = paste0("result2/",cancer_name,'/monocle3_all.RData'))
