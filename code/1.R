# #########首先事要读取癌症的名字
# 
# cancer_names <- list.files("H:/PCD_result/monocle3_r/")
# #########癌症的地址
# getwd()
# cancer_file  <- list.files("H:/PCD_result/结果一/")
# cancer_file <- paste0("H:/PCD_result/结果一/",
#                       cancer_file,"/",cancer_file)
# setwd("D:/R-daoma/pcd_databsee/")
# #########把程序性死亡基因 整理
# Pcd_gene <- read.csv("死亡基因文件.csv")
# Pcd_gene <- Pcd_gene[,2:16]
# colnames(Pcd_gene)
# ############# 单个死亡基因提取
# Pyroptosis <- unique(Pcd_gene$Pyroptosis)
# Pyroptosis <- Pyroptosis[-which(Pyroptosis=="")]
# 
# Ferroptosis <- unique(Pcd_gene$Ferroptosis)
# Ferroptosis <- Ferroptosis[-which(Ferroptosis=="")]
# 
# Autophagy <- unique(Pcd_gene$Autophagy)
# Autophagy <- Autophagy[-which(Autophagy=="")]
# 
# Necroptosis <- unique(Pcd_gene$Necroptosis)
# Necroptosis <- Necroptosis[-which(Necroptosis=="")]
# 
# Cuproptosis <- unique(Pcd_gene$Cuproptosis)
# Cuproptosis <- Cuproptosis[-which(Cuproptosis=="")]
# 
# Parthanatos <- unique(Pcd_gene$Parthanatos)
# Parthanatos <- Parthanatos[-which(Parthanatos=="")]
# 
# Entotie_cell_death <- unique(Pcd_gene$Entotie.cell.death)
# Entotie_cell_death <- Entotie_cell_death[-which(Entotie_cell_death=="")]
# 
# Netotic_cell_death <- unique(Pcd_gene$Netotic.cell.death)
# Netotic_cell_death <- Netotic_cell_death[-which(Netotic_cell_death=="")]
# 
# Lysosome_dependent_cell_death <- unique(Pcd_gene$Lysosome.dependent.cell.death)
# Lysosome_dependent_cell_death <- Lysosome_dependent_cell_death[-which(Lysosome_dependent_cell_death=="")]
# 
# Alkaliptosis <-  unique(Pcd_gene$Alkaliptosis)
# Alkaliptosis <- Alkaliptosis[-which(Alkaliptosis=="")]
# 
# Oxeiptosis  <- unique(Pcd_gene$Oxeiptosis)
# Oxeiptosis <- Oxeiptosis[-which(Oxeiptosis=="")]
# 
# Immunogenic_cell_death <- unique(Pcd_gene$Immunogenic_cell_death)
# Immunogenic_cell_death <- Immunogenic_cell_death[-which(Immunogenic_cell_death=="")]
# 
# Intrinsic_apoptosis <-  unique(Pcd_gene$intrinsic.apoptosis)
# Intrinsic_apoptosis <- Intrinsic_apoptosis[-which(Intrinsic_apoptosis=="")]
# 
# Extrinsic_apoptosis <-  unique(Pcd_gene$extrinsic.apoptosis)
# Extrinsic_apoptosis <- Extrinsic_apoptosis[-which(Extrinsic_apoptosis=="")]
# 
# Anoikis <-  unique(Pcd_gene$Anoikis)
# Anoikis <- Anoikis[-which(Anoikis=="")]
# 
# marker_gene_list <- list(Pyroptosis=Pyroptosis,
#                          Ferroptosis=Ferroptosis,
#                          Autophagy=Autophagy,Necroptosis=Necroptosis,Cuproptosis=Cuproptosis,Parthanatos=Parthanatos,
#                          Entotie_cell_death=Entotie_cell_death,Netotic_cell_death=Netotic_cell_death,Lysosome_dependent_cell_death=Lysosome_dependent_cell_death,
#                          Alkaliptosis=Alkaliptosis,Oxeiptosis=Oxeiptosis,Immunogenic_cell_death=Immunogenic_cell_death,Intrinsic_apoptosis=Intrinsic_apoptosis
#                          ,Extrinsic_apoptosis=Extrinsic_apoptosis,Anoikis=Anoikis)
# 
# 
# ls()
# save(Alkaliptosis,Autophagy,Anoikis,Cuproptosis,
#      Entotie_cell_death,Extrinsic_apoptosis,Ferroptosis,
#      Immunogenic_cell_death,Intrinsic_apoptosis,Lysosome_dependent_cell_death,
#      Necroptosis,Netotic_cell_death,Oxeiptosis,Parthanatos,Pyroptosis,file = "PCD_gene_set.RData")


######################################################################################################
rm(list=ls())
#############读取PCD 
cancer_name <- "UVM"
dir.create(paste0("result1/",cancer_name))
load("PCD_gene_set.RData")
marker_gene_list <- list(Pyroptosis=Pyroptosis,
                         Ferroptosis=Ferroptosis,
                         Autophagy=Autophagy,Necroptosis=Necroptosis,Cuproptosis=Cuproptosis,Parthanatos=Parthanatos,
                         Entotie_cell_death=Entotie_cell_death,Netotic_cell_death=Netotic_cell_death,Lysosome_dependent_cell_death=Lysosome_dependent_cell_death,
                         Alkaliptosis=Alkaliptosis,Oxeiptosis=Oxeiptosis,Immunogenic_cell_death=Immunogenic_cell_death,Intrinsic_apoptosis=Intrinsic_apoptosis
                         ,Extrinsic_apoptosis=Extrinsic_apoptosis,Anoikis=Anoikis)

PCD_names <- c("Alkaliptosis","Anoikis","Autophagy","Cuproptosis",
"Entotie_cell_death","Extrinsic_apoptosis","Ferroptosis",
"Immunogenic_cell_death","Intrinsic_apoptosis",'Lysosome_dependent_cell_death',
"Necroptosis","Netotic_cell_death","Oxeiptosis","Parthanatos","Pyroptosis")
############读取癌症
library(Seurat)
library(tidyverse)
library(GSVA)
library(plyr)
library(ggplot2)
library(tidydr)
library(dplyr)
library(cols4all)
library(edgeR)
library(DESeq2)
library(fgsea)
library(enrichplot)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(magick)
############load 癌症数据
load("H:/PCD_result/结果一/ESCC_data.RData/ESCC_data.RData")
############对癌症细胞的Idents 进行赋值
if(dim(new_data)[2]>10000){#########大于10000细胞固定0.8。小于10000固定0.5
  Idents(new_data) <- "RNA_snn_res.0.8"###########这里需要注意，整合的数据叫做 integrated_snn_res.0.8
  temp <- new_data@meta.data$RNA_snn_res.0.8##############integrated_snn_res.0.8
  new_data <-  AddMetaData(new_data,temp,col.name = "NEW.meta")
}else{
  Idents(new_data) <- "RNA_snn_res.0.5"#######integrated_snn_res.0.5
  temp <- new_data@meta.data$RNA_snn_res.0.5###########integrated_snn_res.0.5
  new_data <-  AddMetaData(new_data,temp,col.name = "NEW.meta")
  }

############寻找all marker gene 
data_conserves_marker <- FindAllMarkers(new_data,
                                        logfc.threshold = 0.5)
########### 超几何分析
the_enrich_chaojihe_ALL <- function(data = data,marker = marker,cluster_number){
  p1 <- c()
  aaaa <- which(colnames(marker) == "gene")
  print(aaaa)
  for(i in 0:cluster_number){
    print(i)
    print(length(intersect(data,marker[which(marker$cluster==i),aaaa])))
    temp <- 1 - phyper(length(intersect(data,marker[which(marker$cluster==i),aaaa])),
                       length(data),
                       22000- length(data),
                       length(marker[which(marker$cluster==i),aaaa]))
    #print(temp)
    p1 <- c(p1, temp)
  }
  p1 <- p.adjust(p1,method = p.adjust.methods,n=length(p1))
  return(p1)
}



#table(new_data@meta.data$RNA_snn_res.0.8)

cluster_number <- max(as.numeric(levels(Idents(new_data))))

#############计算每个PCD的富集得分


chaojihe_enrich_p <- list()
for(i in PCD_names){
  temp <- the_enrich_chaojihe_ALL(get(i),marker = data_conserves_marker,
                                  cluster_number = cluster_number)
  chaojihe_enrich_p[[i]] <- temp}

chaojihe_enrich_p <- do.call(cbind,chaojihe_enrich_p)
rownames(chaojihe_enrich_p) <- 0:cluster_number
gc()
############路径
file_chaojihe <-  paste("result1",cancer_name,"chaojihe.csv",sep = "/")
write.csv(chaojihe_enrich_p,file =file_chaojihe)

#############进行ssGSEA的富集分分析
subpathway <- list(
  Alkaliptosis=Alkaliptosis,Anoikis=Anoikis,
  Autophagy=Autophagy,
  Cuproptosis=Cuproptosis,
  Entotie_cell_death=Entotie_cell_death,
  Extrinsic_apoptosis= Extrinsic_apoptosis,
  Ferroptosis=Ferroptosis,
  Immunogenic_cell_death=Immunogenic_cell_death,
  Intrinsic_apoptosis=Intrinsic_apoptosis,
  Lysosome_dependent_cell_death=Lysosome_dependent_cell_death,
  Necroptosis=Necroptosis,Netotic_cell_death=Netotic_cell_death,
  Oxeiptosis=Oxeiptosis,Parthanatos=Parthanatos,Pyroptosis=Pyroptosis)

#########运行ssGSEA的分析
gc()
ES_list <- list()
for (i in 0:cluster_number){
  print(i)
  cluster_one <- subset(new_data, NEW.meta == i)
  dddd <- which(colnames(data_conserves_marker)=="gene")
  cluster_ono_exp <- as.data.frame(cluster_one@assays$RNA@data)[data_conserves_marker[which(data_conserves_marker$cluster==i),dddd],]
  #计算全部子通路的ssGSEA打分
  if(dim(cluster_ono_exp)[1] <10){
    ES_list[[i+1]] <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)%>% t() %>% as.data.frame() 
    names( ES_list[[i+1]]) <- PCD_names
    next
  }
  

  GBM_ssGSEA_result_all <- gsva(as.matrix(unique(cluster_ono_exp)), subpathway,method='ssgsea',
                                kcdf='Gaussian',abs.ranking=F,
                                min.sz = 1,parallel=8, mx.diff=TRUE)
  ES_list[[i+1]] <- as.data.frame(t(apply(GBM_ssGSEA_result_all,1,mean)))
}

gc()
aaaaaa <-do.call(rbind.fill,ES_list)
#############判断缺失数量
if(dim(aaaaaa)[2]<15){
  PCD_names_bu <-  PCD_names[ !PCD_names%in%  colnames(aaaaaa)]
  for(i in PCD_names_bu){
    aaaaaa[,i] <- NA
  }
}
aaaaaa <-aaaaaa[,PCD_names]
aaaaaa[is.na(aaaaaa)] <- 0
rownames(aaaaaa) <- 0:cluster_number
enrich_file <-  paste("result1",cancer_name,"GSEA_enrich.csv",sep = "/")

write.csv(aaaaaa,file = enrich_file)

######  亚群重注释：
the_ident_1 <- c()
for(i in 0:cluster_number){
  print(i)
  temp <- which(chaojihe_enrich_p[i+1,]<0.05)
  if(length(temp) >1){
    temp_1 <- abs(aaaaaa[i+1,temp])
    
    temp_2 <- max(temp_1[1,])
    temp_3 <- colnames(temp_1)[which(temp_1[1,]==temp_2)]
    if(length(temp_3)>1){
      temp_3 <- temp_3[order(temp_3,decreasing = F)]
      temp_3 <- paste(temp_3,collapse = ",")
      if(temp_3=="Extrinsic_apoptosis,Intrinsic_apoptosis"){
        temp_3 <- "Extrinsic_Intrinsic_apoptosis"
      }
    }
    the_ident_1 <- c(the_ident_1,temp_3) 
  }else if(length(temp)  == 1){
    temp_3 <- names(temp)
    the_ident_1 <- c(the_ident_1,temp_3)
  }else if(length(temp)  == 0){
    the_ident_1 <- c(the_ident_1,"None")
  }
}
names(the_ident_1) <- 0:cluster_number
new_data <- RenameIdents(new_data, the_ident_1)
save(new_data,file =paste("result1",cancer_name,"PCD_ANN.RData",sep = "/") )
###########

###########整体pcd的颜色


# color_PCD <- c("Alkaliptosis"="#f1ca17","Autophagy"="#77CEF3",
#                "Cuproptosis"="#59BD7E",
# "Entotie_cell_death"="#F8C0CB","Extrinsic_apoptosis"="#11659a",
# "Ferroptosis"="#f8c387",
# "Immunogenic_cell_death"="#CA8E4F","Intrinsic_apoptosis"="#B56B52",
# 'Lysosome_dependent_cell_death'="#c06f98",
# "Necroptosis"="#EE856D","Netotic_cell_death"="#b9dec9",
# "Oxeiptosis"="#3BC1C5",
# "Parthanatos"="#813c85","Pyroptosis" = "#f07c82",
# "Extrinsic_apoptosis,Intrinsic_apoptosis"=="#EE2C2C",
# "None"="#8A959B")
# 
# # 
#  "#A52A2A","#41AB5D","#8DD3C7","#B3DE69","#BC80BD","#BEBADA","#C6DBEF","#80B1D3","#FB8072","#EF9BB5","#F3B161"
# # 
# "#539A3F","#EB9B98","#ABCADE","#644092","#B7DA90","#3B77AC","#A35D34","#ED8433","#F1BE7B","#C3B0D0","#CF352B"


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


order_ann <-names( color_PCD)[-17]


###### Plot tsne
tsne_file  <-  paste("result1",cancer_name,"tsne_plot.pdf",sep = "/")
pdf(file = tsne_file,width = 11.69,height = 8.27)
Idents(new_data)
DimPlot(new_data, label = F,reduction = "tsne",repel = T,
        order =rev(c(order_ann,"None")),cols = color_PCD,
        pt.size = 1.5 )+
  theme(legend.position = "right",axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        plot.title = element_blank())+
  theme_dr(xlength = 0.2, #x轴长度
           ylength= 0.2, #y轴长度
           arrow= grid::arrow(length = unit(0.1,"inches"), #箭头大小/长度
                              ends= 'last', type = "closed")) + #箭头描述信息
  theme(panel.grid = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 4.5)))
dev.off()
gc()
########### plot UMAP
umap_file  <-  paste("result1",cancer_name,"umap_plot.pdf",sep = "/")
pdf(file = umap_file,width = 11.69,height = 8.27)
DimPlot(new_data, label = F,reduction = "umap",repel = T,
        order  =rev(c(order_ann,"None")),cols = color_PCD,
        pt.size = 1.5 )+
  theme(legend.position = "right",axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        plot.title = element_blank())+
  theme_dr(xlength = 0.2, #x轴长度
           ylength= 0.2, #y轴长度
           arrow= grid::arrow(length = unit(0.1,"inches"), #箭头大小/长度
                              ends= 'last', type = "closed")) + #箭头描述信息
  theme(panel.grid = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 4.5)))
dev.off()

#######各个pcd tu '#F68282' "grey"

for(i in unique(Idents(new_data))){
  col <- c("#F68282","grey")
  names(col)[1] <- i
 p1 <-  DimPlot(new_data, label = F,reduction = "umap",repel = T,
                order  =rev(c(order_ann,"None")),cols = col,
                pt.size = 1.5 )+
    theme(legend.position = "right",axis.title.y=element_blank(), 
          axis.text.y=element_blank(),
          axis.title.x=element_blank(), 
          axis.text.x=element_blank(),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          plot.title = element_blank())+
    theme_dr(xlength = 0.2, #x轴长度
             ylength= 0.2, #y轴长度
             arrow= grid::arrow(length = unit(0.1,"inches"), #箭头大小/长度
                                ends= 'last', type = "closed")) + #箭头描述信息
    theme(panel.grid = element_blank())+
    guides(color = guide_legend(override.aes = list(size = 4.5)))
 p2 <-  DimPlot(new_data, label = F,reduction = "tsne",repel = T,
                order =rev(c(order_ann,"None")),cols = col,
                pt.size = 1.5 )+
   theme(legend.position = "right",axis.title.y=element_blank(), 
         axis.text.y=element_blank(),
         axis.title.x=element_blank(), 
         axis.text.x=element_blank(),
         axis.line=element_blank(),
         axis.ticks=element_blank(),
         plot.title = element_blank())+
   theme_dr(xlength = 0.2, #x轴长度
            ylength= 0.2, #y轴长度
            arrow= grid::arrow(length = unit(0.1,"inches"), #箭头大小/长度
                               ends= 'last', type = "closed")) + #箭头描述信息
   theme(panel.grid = element_blank())+
   guides(color = guide_legend(override.aes = list(size = 4.5)))
 
 
 p1_file  <-  paste0("result1/",cancer_name,"/",i,"_umap_plot.pdf")
 pdf(file = p1_file,width = 11.69,height = 8.27)
 print(p1)
 dev.off()
 
 p2_file  <-  paste0("result1/",cancer_name,"/",i,"_tsne_plot.pdf")
 pdf(file = p2_file,width = 11.69,height = 8.27)
 print(p2)
 dev.off()
}

#######
####### 每个癌症，每个死亡的Hallmarker 富集分析；
new_data <- AddMetaData(new_data,metadata = Idents(new_data),col.name = "active.ident")
Cell <- unique(new_data$active.ident)
group_data <- new_data@meta.data
group_data <- data.frame(ID = rownames(group_data),
                         Group = group_data$active.ident) 
difF_list <- list()
temp_data <- Cell
data_GSEA <- data.frame(Description=0,enrichmentScore=0,NES=0,p.adjust=0,Cell=0 ,cancer=0)
gc()
for(z in 1:length(Cell)){
  print(Cell[z])
  res <- FindMarkers(new_data,ident.1 = Cell[z])
  alldiff <- res%>%na.omit()
  alldiff$type <- ifelse(alldiff$p_val_adj>0.05,'No-Sig',
                         ifelse(alldiff$avg_log2FC>1,'Up',
                                ifelse(alldiff$avg_log2FC< -1,'Down','No-Sig')))
  #接下来开始重头戏，开始画GSEA table
  ## 根据logfc降序排列基因
  alldiff <- alldiff[order(alldiff$avg_log2FC,decreasing = T),]
  ## fgsea中输入的关键基因信息
  id <- alldiff$avg_log2FC
  names(id) <- rownames(alldiff)
  gmtfile <- "./h.all.v7.5.1.symbols.gmt"
  hallmark <- read.gmt(gmtfile)
  hallmark$term <- gsub('HALLMARK_','',hallmark$term)
  hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)
  fgseaRes <- fgsea(pathways = hallmark.list, 
                    stats = id,
                    minSize=1,
                    maxSize=10000,
                    nperm=1000)
  sig <- fgseaRes[fgseaRes$padj<0.05,]
  sig <- sig[order(sig$NES,decreasing = T),]
  
  difF_list[[z]] <- sig
  gc()
}

for(f in 1:length(difF_list)){
  temp_11 <- as.data.frame(difF_list[[f]])
  if(dim(temp_11)[1]==0){next()}
  temp_22 <- data.frame(Description= temp_11$pathway, enrichmentScore = temp_11$ES,
                        NES = temp_11$NES,  p.adjust = temp_11$padj
                        , Cell = Cell[f],cancer = cancer_name)
  
  
  data_GSEA <- rbind(data_GSEA,temp_22)
}
data_GSEA <- data_GSEA[-1,]
write.csv(data_GSEA,file = paste0("result1/",cancer_name,"/","GSEA_HM.csv"))

#########画KEGG通路图：
load("H:/PCD_database/result1/UVM/ALL.RData")
cell_name_states <- data.frame(cell_names = colnames(new_data),
                               states = Idents(new_data))
######首先得到这个癌症有哪些程序性死亡状态
PCD_names_this_cancer <-  unique(Idents(new_data))
PCD_names_this_cancer <- PCD_names_this_cancer[PCD_names_this_cancer %in% c("Intrinsic_apoptosis",
                                                                            "Extrinsic_apoptosis",
                                                                            "Extrinsic_Intrinsic_apoptosis",
                                                                            "Autophagy","Ferroptosis",
                                                                            "Pyroptosis","Necroptosis")]
PCD_names_this_cancer <- as.character(PCD_names_this_cancer)
######根据程序性死亡的种类进行绘图 5种
###### Apoptosis 的细胞进行加和平均
#提取这一种类的细胞名字
getwd()
setwd(paste0("H:/PCD_database/result1/",cancer_name))########修改一下路径
temp <- length(PCD_names_this_cancer)

marker_gene_list <- list(Pyroptosis=Pyroptosis,
                         Ferroptosis=Ferroptosis,
                         Autophagy=Autophagy,Necroptosis=Necroptosis,Cuproptosis=Cuproptosis,Parthanatos=Parthanatos,
                         Entotie_cell_death=Entotie_cell_death,Netotic_cell_death=Netotic_cell_death,Lysosome_dependent_cell_death=Lysosome_dependent_cell_death,
                         Alkaliptosis=Alkaliptosis,Oxeiptosis=Oxeiptosis,Immunogenic_cell_death=Immunogenic_cell_death,Intrinsic_apoptosis=Intrinsic_apoptosis
                         ,Extrinsic_apoptosis=Extrinsic_apoptosis,Anoikis=Anoikis,
                         Extrinsic_Intrinsic_apoptosis=c(unique(c(Extrinsic_apoptosis,Intrinsic_apoptosis))))




sort(names(marker_gene_list))
gc()
for(d in 1:temp){
  print(d)
  ########apoptosis 
  cells_temp_names <- cell_name_states %>% filter(states == PCD_names_this_cancer[d])
  #提取所有的基因
  gene_temp <-marker_gene_list[[PCD_names_this_cancer[d]]][marker_gene_list[[PCD_names_this_cancer[d]]] %in% 
                                                 rownames(new_data)] 
  ##加和计算，并同时进行数据的整理
  id_list <- mapIds(org.Hs.eg.db,gene_temp,"ENTREZID","SYMBOL")
  data_caner_temp <- as.data.frame(new_data@assays$RNA@data)
  #ddd <- as.data.frame(new_data@assays$RNA@scale.data)
  
  data_caner_temp <- data_caner_temp[gene_temp,cells_temp_names$cell_names]
  gene_expression_mean <- rowMeans(as.matrix(data_caner_temp),
                                   na.rm = T) 
  gene_expression_mean <-as.data.frame(gene_expression_mean)
  gene_expression_mean <- gene_expression_mean[names(id_list),]
  gene_expression_mean <- data.frame(ID_entrezid = id_list,gene_expression_mean )
  gene_expression_mean <- na.omit(gene_expression_mean)
  rownames(gene_expression_mean) <- gene_expression_mean$ID_entrezid
  gene_expression_mean <- as.matrix(gene_expression_mean)
  gene_expression_mean[,2] <- as.numeric(gene_expression_mean[,2])
  gene_expression_mean_1 <- as.numeric(gene_expression_mean[,2])
  names(gene_expression_mean_1) <- gene_expression_mean[,1]
  if(PCD_names_this_cancer[d] == "Intrinsic_apoptosis" | PCD_names_this_cancer[d] == "Extrinsic_apoptosis" |
     PCD_names_this_cancer[d] == "Extrinsic_Intrinsic_apoptosis"){
    p <- pathview(gene.data =gene_expression_mean_1, pathway.id = "04210", species = "hsa",
                  out.suffix = paste0(cancer_name,"_",PCD_names_this_cancer[d],"_KEGG"), kegg.native = T, 
                  same.layer = F,min.nnodes = 1,res=350)
  }else if(PCD_names_this_cancer[d] == "Autophagy"){
    p <- pathview(gene.data =gene_expression_mean_1, pathway.id = "04140", species = "hsa",
                  out.suffix = paste0(cancer_name,"_",PCD_names_this_cancer[d],"_KEGG"), kegg.native = T, 
                  same.layer = F,min.nnodes = 1,res=350)
  }else if(PCD_names_this_cancer[d] == "Ferroptosis"){
    p <- pathview(gene.data =gene_expression_mean_1, pathway.id = "04216", species = "hsa",
                  out.suffix = paste0(cancer_name,"_",PCD_names_this_cancer[d],"_KEGG"), kegg.native = T, 
                  same.layer = F,min.nnodes = 1,res=350)
  }else if(PCD_names_this_cancer[d] == "Pyroptosis"){
    p <- pathview(gene.data =gene_expression_mean_1, pathway.id = "04621", species = "hsa",
                  out.suffix = paste0(cancer_name,"_",PCD_names_this_cancer[d],"_KEGG"), kegg.native = T, 
                  same.layer = F,min.nnodes = 1,res=350)
  }else if(PCD_names_this_cancer[d] == "Necroptosis"){
    p <- pathview(gene.data =gene_expression_mean_1, pathway.id = "04217", species = "hsa",
                  out.suffix = paste0(cancer_name,"_",PCD_names_this_cancer[d],"_KEGG"), kegg.native = T, 
                  same.layer = F,min.nnodes = 1,res=350)
  }
}

setwd("../")
setwd("../")

gc()

save.image(file = paste0("D:\\R-daoma\\pcd_databsee\\result1\\",cancer_name,"\\ALL.RData"))#####修改路径
