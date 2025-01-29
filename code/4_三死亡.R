rm(list=ls())
library(readxl)
library(readr)
library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(corrplot)
library(data.table)
library(igraph)
library(monocle3)
#############结果四的分析代码
load("互作统计.RData")##############需要看一下
############################一对重要互作的癌症的代码-----------------
#############对癌症进行一个KEGG和通路网络构建
cancer_name <- "ALL"
#dir.create("H:/PCD_database/result4")
dir.create(paste0("H:/PCD_database/result4/",cancer_name))
interaction_data <- as.data.frame(interaction_data)
temp <-unique(c( as.character(interaction_data[which(interaction_data$cancer==cancer_name),3]),
                 as.character(interaction_data[which(interaction_data$cancer==cancer_name),4]) ))

A_death_name <- temp[1]
B_death_name <- temp[2]
C_death_name <- temp[3]
###########这个部分需要联合结果一，结果二，结果三的结果和数据=
###########首先是读入结果一
data_file <- paste0("H:/PCD_database/result1/",cancer_name,"/PCD_ANN.RData")
load(data_file)

###########对结果一数据进行差异基因的一个计算
DEG_new <- FindAllMarkers(new_data,logfc.threshold = 0.25)
############提取一下重要亚群的差异基因
DEG_new <- DEG_new$gene[DEG_new$cluster %in%temp ]
DEG_new <- unique(DEG_new)
rm(new_data)

############提取结果二中的基因数据-----------
data_file <- paste0("H:/PCD_database/result2/",cancer_name,"/gene_module.RData")
load(data_file)
data_file <- paste0("H:/PCD_database/result2/",cancer_name,"/CDS.RData")
load(data_file)
############对模块数据进行一个整理
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$Death_state)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
############ A death 高表达的module
A_module <-rownames(agg_mat)[order(agg_mat[,A_death_name],decreasing = T)[1:3]]
############ B death 高表达的module
B_module <-rownames(agg_mat)[order(agg_mat[,B_death_name],decreasing = T)[1:3]]
########### C death 高表达的module
C_module <-rownames(agg_mat)[order(agg_mat[,C_death_name],decreasing = T)[1:3]]
rm(cds)
############提取A B modules 的基因
module_number <- as.numeric(str_split_fixed(c(A_module,B_module,C_module)," ",2)[,2])
module_gene <- gene_module_df$id[gene_module_df$module %in% module_number]
module_gene <- unique(module_gene)

#############提取结果三的数据
############# 在这一部 我们使用significant 的means
#############读取数据和简单处理
temp <- grep("statistical_analysis_significant_means",
             list.files(paste0("H:/PCD_database/result3_database/",
                               cancer_name,"/result_cpdb")),
             value = T)
data_file <- paste0("H:/PCD_database/result3_database/",cancer_name,"/result_cpdb/",temp)

interaction_means <- read.table(data_file,
                                header = T,fill = T,
                                sep = "\t")

colnames(interaction_means) <- gsub(pattern = "\\.",replacement = "|",colnames(interaction_means))
colnames(interaction_means) <- gsub(pattern = "nan",replacement = "None",colnames(interaction_means))
colnames(interaction_means)
interaction_means <- interaction_means %>% filter(directionality=="Ligand-Receptor")
#############首先需要写出我们的互作
interaction_A_B <- paste0(A_death_name,"|",B_death_name)
interaction_B_A <- paste0(B_death_name,"|",A_death_name)
interaction_B_C <- paste0(B_death_name,"|",C_death_name)
interaction_C_B <- paste0(C_death_name,"|",B_death_name)
interaction_A_C <- paste0(A_death_name,"|",C_death_name)
interaction_C_A <- paste0(C_death_name,"|",A_death_name)

number <- unique(c(which(interaction_means[,interaction_A_B] >0),
                   which(interaction_means[,interaction_B_A] >0),
                   which(interaction_means[,interaction_B_C] >0),
                   which(interaction_means[,interaction_C_B] >0),
                   which(interaction_means[,interaction_A_C] >0),
                   which(interaction_means[,interaction_C_A] >0)))

ligand_receptor_gene <- unique(c(interaction_means$gene_a[number],
                                 interaction_means$gene_b[number]))
ligand_receptor_gene <-ligand_receptor_gene[-which(ligand_receptor_gene=="")] 

#################最后综合所得到的三个基因集合
gene_set <- c(DEG_new,module_gene,ligand_receptor_gene) 
gene_set <- unique(gene_set)
save(gene_set,file = paste0("H:/PCD_database/result4/",cancer_name,"/gene_set.RData"))
#load("H:/PCD_database/result4/ALL/gene_set.RData")
#################计算KEGG的一个通路富集
gene_set <- bitr(gene_set,fromType = "SYMBOL",
                 toType = c("ENSEMBL","ENTREZID"),
                 OrgDb = "org.Hs.eg.db")

########富集分析
gene_set_enrich_kegg <- enrichKEGG(gene = gene_set$ENTREZID,
                                   organism = "hsa",
                                   keyType = "kegg",          
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05)

save(gene_set_enrich_kegg ,file = paste0("H:/PCD_database/result4/",
                                         cancer_name,"/gene_set_enrich_kegg.RData"))
#load("H:/PCD_database/result4/ALL/gene_set_enrich_kegg.RData")
##########筛选P值
gene_set_enrich_kegg <-  as.data.frame(gene_set_enrich_kegg@result)

gene_set_enrich_kegg <- gene_set_enrich_kegg %>% filter(p.adjust <0.05)


gene_set_enrich_kegg$Description <- factor(gene_set_enrich_kegg$Description, 
                                           levels = rev(gene_set_enrich_kegg$Description))

if(dim(gene_set_enrich_kegg)[1]>20){
  gene_set_enrich_kegg_plot <- gene_set_enrich_kegg[1:20,]
}else{
  gene_set_enrich_kegg_plot <- gene_set_enrich_kegg
}
#KEGG富集柱形图绘制

#自定义主题：
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))

#创建一个可套用的绘图函数减少重复代码：
KEGG_bar <- function(x){
  y <- get(x)
  ggplot(data = y, aes(x = Count, y = Description, fill = -log10(p.adjust))) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_y_discrete(labels = function(y) 
      str_wrap(y, width = 30) ) + #label换行，部分term描述太长
    labs(x = "Gene Number", y = "",
         title = "KEGG enrichment barplot") +
    theme_bw() +
    mytheme+
    theme(axis.title.x=element_text(color = "black"),
          axis.title.y=element_text(color = "black"),
          axis.text.y = element_text(color = "black",size = 15),
          axis.text.x = element_text(color = "black",size = 15))+
    theme(legend.text = element_text(color = "black"),
          legend.title =element_text(color = "black") )
}

p3 <- KEGG_bar("gene_set_enrich_kegg_plot")+scale_fill_distiller(palette = "Oranges",direction = 1)


pdf_file <- paste0("H:/PCD_database/result4/",cancer_name,"/","KEGG_bar.pdf")
pdf(file = pdf_file,width = 8.27,height = 11.67)
print(p3)
dev.off()

########################################网络计算前期准备
#########读入PPI数据
load("ppi_yunxing.RData")
ppi[,3:4] <- lapply(ppi[,3:4],as.character) #将PPI网络的entrez列转为字符型

#########列出KEGG所有的通路组合
path <- gene_set_enrich_kegg[,c(3,4,10)]
all_path <- NULL
for(i in 1:(nrow(path)-1)){
  for (j in (i+1):nrow(path)){
    xx <- cbind(path[i,c(1:3)],path[j,c(1:3)])
    colnames(xx) <- ""
    all_path <- rbind(all_path,xx)
  }
}
all_path <- all_path[,c(1,4,2,5,3,6)]
colnames(all_path) <- c("id1","id2","name1","name2","gene1","gene2")
rownames(all_path) <- c()

###################确定两个death名字
A_death_name
B_death_name
C_death_name
##################再次读入癌症数据
data_file <- paste0("H:/PCD_database/result1/",cancer_name,"/PCD_ANN.RData")
load(data_file)
PPI <- ppi
####################### A VS B
#################提取细胞名和分组
group_info <- data.frame(cell_names = colnames(new_data),
                         group = Idents(new_data))

group_info <-group_info[ group_info$group %in% c(A_death_name,B_death_name),]
#################提取细胞表达
new_data_AB <- new_data[,group_info$cell_names]
exp <- as.matrix(new_data_AB@assays$RNA@data)
name <- bitr(rownames(exp),fromType = "SYMBOL",toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
exp <- exp[name$SYMBOL,]
rownames(exp) <- name$ENTREZID
rm(new_data_AB)
gc()

#################确定A_death_name和B_death_name的细胞位置
A_death_local <-  which(group_info$group == A_death_name)
B_death_local <-  which(group_info$group == B_death_name)

#########
#########
#计算通路之间的crosstalk
sum3<-c()
for(i in 1:dim(all_path)[1]){
  print(i)
  p1 <- unlist(strsplit(all_path[i,5],"\\/"))
  p2 <- unlist(strsplit(all_path[i,6],"\\/"))
  
  g_id1 <- which(PPI[,3] %in% p1)
  g_id2 <- which(PPI[,4] %in% p2)
  id3 <- intersect(na.omit(g_id1),na.omit(g_id2))
  sum1 <- 0
  if(length(id3)!=0){
    for (j in 1:length(id3)) {
      x <- PPI$ENTREZID1[id3[j]]
      y <- PPI$ENTREZID2[id3[j]]
      if(x%in%rownames(exp) & y%in%rownames(exp)){
        temp1 <- as.data.frame(table(exp[x,A_death_local]==exp[x,B_death_local]))
        temp2 <-  as.data.frame(table(exp[y,A_death_local]==exp[y,B_death_local]))
        if(dim(temp1)[1]<2){
          if(temp1[1,1] == "FALSE"){
            temp1 <- rbind(temp1,
                           data.frame(Var1 ="TRUE",
                                      Freq = 0))
          }else{
            temp1 <- rbind(temp1,
                           data.frame(Var1 ="FALSE",
                                      Freq = 0))
          }
        }
        if(dim(temp2)[1]<2){
          if(temp2[1,1] == "FALSE"){
            temp2 <- rbind(temp2,
                           data.frame(Var1 ="TRUE",
                                      Freq = 0))
          }else{
            temp2 <- rbind(temp2,
                           data.frame(Var1 ="FALSE",
                                      Freq = 0))
          }
        }
        rownames(temp1) <- temp1[,1]
        rownames(temp2) <- temp1[,1]
        
        if((temp1["TRUE",2] == length(B_death_local))|
           (temp2["TRUE",2] == length(B_death_local))){
          s1=0}else{
            s1 <- (-1)*log(t.test(exp[x,A_death_local],exp[x,B_death_local],alternative = "two.sided")$p.value)+
              (-1)*log(t.test(exp[y,A_death_local],exp[y,B_death_local],alternative = "two.sided")$p.value)+
              (-1)*log(cor.test(as.numeric(exp[x,]),as.numeric(exp[y,]),alternative = "two.sided")$p.value)}
      }
      else{s1=0}
      sum1=sum1+s1
    }
  }
  
  g_id3 <- which(PPI[,4] %in% p1)
  g_id4 <- which(PPI[,3] %in% p2)
  id4 <- intersect(na.omit(g_id3),na.omit(g_id4))
  id5 <- setdiff(id4,id3)
  sum2 <- 0
  if(length(id5)!=0){
    for (z in 1:length(id5)) {
      m <- PPI$ENTREZID2[id5[z]]
      n <- PPI$ENTREZID1[id5[z]]
      if(m%in%rownames(exp) & n%in%rownames(exp)){
        temp1 <- as.data.frame(table(exp[m,A_death_local]==exp[m,B_death_local]))
        temp2 <-  as.data.frame(table(exp[n,A_death_local]==exp[n,B_death_local]))
        if(dim(temp1)[1]<2){
          if(temp1[1,1] == "FALSE"){
            temp1 <- rbind(temp1,
                           data.frame(Var1 ="TRUE",
                                      Freq = 0))
          }else{
            temp1 <- rbind(temp1,
                           data.frame(Var1 ="FALSE",
                                      Freq = 0))
          }
        }
        if(dim(temp2)[1]<2){
          if(temp2[1,1] == "FALSE"){
            temp2 <- rbind(temp2,
                           data.frame(Var1 ="TRUE",
                                      Freq = 0))
          }else{
            temp2 <- rbind(temp2,
                           data.frame(Var1 ="FALSE",
                                      Freq = 0))
          }
        }
        rownames(temp1) <- temp1[,1]
        rownames(temp2) <- temp1[,1]
        
        
        if((temp1["TRUE",2] == length(B_death_local))|
           (temp2["TRUE",2] == length(B_death_local))){
          s2=0}else{
            s2 <- (-1)*log(t.test(exp[m,A_death_local],exp[m,B_death_local],alternative = "two.sided")$p.value)+
              (-1)*log(t.test(exp[n,A_death_local],exp[n,B_death_local],alternative = "two.sided")$p.value)+
              (-1)*log(cor.test(as.numeric(exp[m,]),as.numeric(exp[n,]),alternative = "two.sided")$p.value)}
      }
      else{s2=0}
      sum2=sum2+s2
    }
  }
  temp <- sum1+sum2
  sum3 <-c(sum3,temp)
}
##############赋值并得到原始的cross
score1 <- sum3
int_path <- all_path[which(score1>0&score1!=Inf),]
score2 <- score1[which(score1>0&score1!=Inf)]
cross_A_B <- cbind.data.frame(int_path,score2)

############ B VS C

#################提取细胞名和分组
group_info <- data.frame(cell_names = colnames(new_data),
                         group = Idents(new_data))

group_info <-group_info[ group_info$group %in% c(B_death_name,C_death_name),]
#################提取细胞表达
new_data_BC <- new_data[,group_info$cell_names]
exp <- as.matrix(new_data_BC@assays$RNA@data)
name <- bitr(rownames(exp),fromType = "SYMBOL",toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
exp <- exp[name$SYMBOL,]
rownames(exp) <- name$ENTREZID
rm(new_data_BC)
gc()

#################确定A_death_name和B_death_name的细胞位置
B_death_local <-  which(group_info$group == B_death_name)
C_death_local <-  which(group_info$group == C_death_name)

#########
#########
#计算通路之间的crosstalk
sum3<-c()
for(i in 1:dim(all_path)[1]){
  print(i)
  p1 <- unlist(strsplit(all_path[i,5],"\\/"))
  p2 <- unlist(strsplit(all_path[i,6],"\\/"))
  
  g_id1 <- which(PPI[,3] %in% p1)
  g_id2 <- which(PPI[,4] %in% p2)
  id3 <- intersect(na.omit(g_id1),na.omit(g_id2))
  sum1 <- 0
  if(length(id3)!=0){
    for (j in 1:length(id3)) {
      x <- PPI$ENTREZID1[id3[j]]
      y <- PPI$ENTREZID2[id3[j]]
      if(x%in%rownames(exp) & y%in%rownames(exp)){
        temp1 <- as.data.frame(table(exp[x,B_death_local]==exp[x,C_death_local]))
        temp2 <-  as.data.frame(table(exp[y,B_death_local]==exp[y,C_death_local]))
        if(dim(temp1)[1]<2){
          if(temp1[1,1] == "FALSE"){
            temp1 <- rbind(temp1,
                           data.frame(Var1 ="TRUE",
                                      Freq = 0))
          }else{
            temp1 <- rbind(temp1,
                           data.frame(Var1 ="FALSE",
                                      Freq = 0))
          }
        }
        if(dim(temp2)[1]<2){
          if(temp2[1,1] == "FALSE"){
            temp2 <- rbind(temp2,
                           data.frame(Var1 ="TRUE",
                                      Freq = 0))
          }else{
            temp2 <- rbind(temp2,
                           data.frame(Var1 ="FALSE",
                                      Freq = 0))
          }
        }
        rownames(temp1) <- temp1[,1]
        rownames(temp2) <- temp1[,1]
        
        if((temp1["TRUE",2] == length(C_death_local))|
           (temp2["TRUE",2] == length(C_death_local))){
          s1=0}else{
            s1 <- (-1)*log(t.test(exp[x,B_death_local],exp[x,C_death_local],alternative = "two.sided")$p.value)+
              (-1)*log(t.test(exp[y,B_death_local],exp[y,C_death_local],alternative = "two.sided")$p.value)+
              (-1)*log(cor.test(as.numeric(exp[x,]),as.numeric(exp[y,]),alternative = "two.sided")$p.value)}
      }
      else{s1=0}
      sum1=sum1+s1
    }
  }
  
  g_id3 <- which(PPI[,4] %in% p1)
  g_id4 <- which(PPI[,3] %in% p2)
  id4 <- intersect(na.omit(g_id3),na.omit(g_id4))
  id5 <- setdiff(id4,id3)
  sum2 <- 0
  if(length(id5)!=0){
    for (z in 1:length(id5)) {
      m <- PPI$ENTREZID2[id5[z]]
      n <- PPI$ENTREZID1[id5[z]]
      if(m%in%rownames(exp) & n%in%rownames(exp)){
        temp1 <- as.data.frame(table(exp[m,B_death_local]==exp[m,C_death_local]))
        temp2 <-  as.data.frame(table(exp[n,B_death_local]==exp[n,C_death_local]))
        if(dim(temp1)[1]<2){
          if(temp1[1,1] == "FALSE"){
            temp1 <- rbind(temp1,
                           data.frame(Var1 ="TRUE",
                                      Freq = 0))
          }else{
            temp1 <- rbind(temp1,
                           data.frame(Var1 ="FALSE",
                                      Freq = 0))
          }
        }
        if(dim(temp2)[1]<2){
          if(temp2[1,1] == "FALSE"){
            temp2 <- rbind(temp2,
                           data.frame(Var1 ="TRUE",
                                      Freq = 0))
          }else{
            temp2 <- rbind(temp2,
                           data.frame(Var1 ="FALSE",
                                      Freq = 0))
          }
        }
        rownames(temp1) <- temp1[,1]
        rownames(temp2) <- temp1[,1]
        
        
        if((temp1["TRUE",2] == length(C_death_local))|
           (temp2["TRUE",2] == length(C_death_local))){
          s2=0}else{
            s2 <- (-1)*log(t.test(exp[m,B_death_local],exp[m,C_death_local],alternative = "two.sided")$p.value)+
              (-1)*log(t.test(exp[n,B_death_local],exp[n,C_death_local],alternative = "two.sided")$p.value)+
              (-1)*log(cor.test(as.numeric(exp[m,]),as.numeric(exp[n,]),alternative = "two.sided")$p.value)}
      }
      else{s2=0}
      sum2=sum2+s2
    }
  }
  temp <- sum1+sum2
  sum3 <-c(sum3,temp)
}
##############赋值并得到原始的cross
score1 <- sum3
int_path <- all_path[which(score1>0&score1!=Inf),]
score2 <- score1[which(score1>0&score1!=Inf)]
cross_B_C <- cbind.data.frame(int_path,score2)

#################### A VS C
#################提取细胞名和分组
group_info <- data.frame(cell_names = colnames(new_data),
                         group = Idents(new_data))

group_info <-group_info[ group_info$group %in% c(A_death_name,C_death_name),]
#################提取细胞表达
new_data_AC <- new_data[,group_info$cell_names]
exp <- as.matrix(new_data_AC@assays$RNA@data)
name <- bitr(rownames(exp),fromType = "SYMBOL",toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
exp <- exp[name$SYMBOL,]
rownames(exp) <- name$ENTREZID
rm(new_data_AC)
gc()

#################确定A_death_name和B_death_name的细胞位置
A_death_local <-  which(group_info$group == A_death_name)
C_death_local <-  which(group_info$group == C_death_name)

#########
#########
#计算通路之间的crosstalk
sum3<-c()
for(i in 1:dim(all_path)[1]){
  print(i)
  p1 <- unlist(strsplit(all_path[i,5],"\\/"))
  p2 <- unlist(strsplit(all_path[i,6],"\\/"))
  
  g_id1 <- which(PPI[,3] %in% p1)
  g_id2 <- which(PPI[,4] %in% p2)
  id3 <- intersect(na.omit(g_id1),na.omit(g_id2))
  sum1 <- 0
  if(length(id3)!=0){
    for (j in 1:length(id3)) {
      x <- PPI$ENTREZID1[id3[j]]
      y <- PPI$ENTREZID2[id3[j]]
      if(x%in%rownames(exp) & y%in%rownames(exp)){
        temp1 <- as.data.frame(table(exp[x,A_death_local]==exp[x,C_death_local]))
        temp2 <-  as.data.frame(table(exp[y,A_death_local]==exp[y,C_death_local]))
        if(dim(temp1)[1]<2){
          if(temp1[1,1] == "FALSE"){
            temp1 <- rbind(temp1,
                           data.frame(Var1 ="TRUE",
                                      Freq = 0))
          }else{
            temp1 <- rbind(temp1,
                           data.frame(Var1 ="FALSE",
                                      Freq = 0))
          }
        }
        if(dim(temp2)[1]<2){
          if(temp2[1,1] == "FALSE"){
            temp2 <- rbind(temp2,
                           data.frame(Var1 ="TRUE",
                                      Freq = 0))
          }else{
            temp2 <- rbind(temp2,
                           data.frame(Var1 ="FALSE",
                                      Freq = 0))
          }
        }
        rownames(temp1) <- temp1[,1]
        rownames(temp2) <- temp1[,1]
        
        if((temp1["TRUE",2] == length(C_death_local))|
           (temp2["TRUE",2] == length(C_death_local))){
          s1=0}else{
            s1 <- (-1)*log(t.test(exp[x,A_death_local],exp[x,C_death_local],alternative = "two.sided")$p.value)+
              (-1)*log(t.test(exp[y,A_death_local],exp[y,C_death_local],alternative = "two.sided")$p.value)+
              (-1)*log(cor.test(as.numeric(exp[x,]),as.numeric(exp[y,]),alternative = "two.sided")$p.value)}
      }
      else{s1=0}
      sum1=sum1+s1
    }
  }
  
  g_id3 <- which(PPI[,4] %in% p1)
  g_id4 <- which(PPI[,3] %in% p2)
  id4 <- intersect(na.omit(g_id3),na.omit(g_id4))
  id5 <- setdiff(id4,id3)
  sum2 <- 0
  if(length(id5)!=0){
    for (z in 1:length(id5)) {
      m <- PPI$ENTREZID2[id5[z]]
      n <- PPI$ENTREZID1[id5[z]]
      if(m%in%rownames(exp) & n%in%rownames(exp)){
        temp1 <- as.data.frame(table(exp[m,A_death_local]==exp[m,C_death_local]))
        temp2 <-  as.data.frame(table(exp[n,A_death_local]==exp[n,C_death_local]))
        if(dim(temp1)[1]<2){
          if(temp1[1,1] == "FALSE"){
            temp1 <- rbind(temp1,
                           data.frame(Var1 ="TRUE",
                                      Freq = 0))
          }else{
            temp1 <- rbind(temp1,
                           data.frame(Var1 ="FALSE",
                                      Freq = 0))
          }
        }
        if(dim(temp2)[1]<2){
          if(temp2[1,1] == "FALSE"){
            temp2 <- rbind(temp2,
                           data.frame(Var1 ="TRUE",
                                      Freq = 0))
          }else{
            temp2 <- rbind(temp2,
                           data.frame(Var1 ="FALSE",
                                      Freq = 0))
          }
        }
        rownames(temp1) <- temp1[,1]
        rownames(temp2) <- temp1[,1]
        
        
        if((temp1["TRUE",2] == length(C_death_local))|
           (temp2["TRUE",2] == length(C_death_local))){
          s2=0}else{
            s2 <- (-1)*log(t.test(exp[m,A_death_local],exp[m,C_death_local],alternative = "two.sided")$p.value)+
              (-1)*log(t.test(exp[n,A_death_local],exp[n,C_death_local],alternative = "two.sided")$p.value)+
              (-1)*log(cor.test(as.numeric(exp[m,]),as.numeric(exp[n,]),alternative = "two.sided")$p.value)}
      }
      else{s2=0}
      sum2=sum2+s2
    }
  }
  temp <- sum1+sum2
  sum3 <-c(sum3,temp)
}

##############赋值并得到原始的cross
score1 <- sum3
int_path <- all_path[which(score1>0&score1!=Inf),]
score2 <- score1[which(score1>0&score1!=Inf)]
cross_A_C <- cbind.data.frame(int_path,score2)

############## 对三个数据进行一个提取和合并######1.0 取交集，希望2.0 取并集
the_rownames_cross <- Reduce(intersect,list(rownames(cross_A_B),
                                             rownames(cross_A_C),
                                             rownames(cross_B_C)))

cross_A_B_1 <- cross_A_B[the_rownames_cross,]
cross_A_C_1 <- cross_A_C[the_rownames_cross,]
cross_B_C_1 <- cross_B_C[the_rownames_cross,]
############计算
cross_A_B_1_score <- cross_A_B_1$score2*0.4
cross_A_C_1_score <- cross_A_C_1$score2*0.2
cross_B_C_1_score <- cross_B_C_1$score2*0.4

########### 得到
cross_final <- all_path[the_rownames_cross,]
cross_final$score2 <- cross_A_B_1_score+cross_A_C_1_score+cross_B_C_1_score
##########写出
data_file <- paste0("H:/PCD_database/result4/",cancer_name,"/cross_raw.txt")
write.table(cross_final,data_file,quote = F, sep = "\t",row.names = F)

#aaaa <- read.table( paste0("H:/PCD_database/result4/",cancer_name,"/cross_raw.txt"))
#############进行cross的筛选和画图
crosstalk_1 <- cross_final
crosstalk_2 <- data.frame(crosstalk_1$name1,crosstalk_1$name2,crosstalk_1$score2)
colnames(crosstalk_2) <- c("Kegg1","Kegg2","score")
crosstalk_3 <- data.frame(crosstalk_1$name2,crosstalk_1$name1,crosstalk_1$score2)
colnames(crosstalk_3) <- c("Kegg1","Kegg2","score")
crosstalk_4 <- rbind(crosstalk_2,crosstalk_3)
aaa <- dcast(crosstalk_4,Kegg1~Kegg2,value.var = "score")
aaa[is.na(aaa)] <- 0
rownames(aaa) <- aaa[,1]
aaa <- aaa[,-1]
aaa<- as.matrix(aaa)
pdf_file <- paste0("H:/PCD_database/result4/",cancer_name,"/","crosstalk_socre.pdf")
pdf(file = pdf_file,width = 25,height = 25)
corrplot(aaa,is.corr = FALSE, 
         type="upper", 
         tl.col="black", tl.srt=50, #修改字体, #显著性筛选
         diag=T,tl.pos = 'td', tl.cex = 0.8, order = "AOE",
         col =rev(COL2("PuOr")),outline = F, cl.pos = 'r')

dev.off()

########################筛选
# if()

########################提取node的名字
########################首先是使用igraph分析
cancer_name <- "ALL"
cross_final <- read.table("H:/PCD_database/result4/ALL/cross_raw.txt",header = T,fill = T,
                       sep = "\t")
load("H:/PCD_database/result4/ALL/gene_set_enrich_kegg.RData")
the_cross_1 <- cross_final[,c(3,4)]
colnames(the_cross_1) <- c("from","to")

g <- graph_from_data_frame(the_cross_1, directed=F)

temp <- names(degree(g,mode = "total")[order(degree(g,mode = "total"),
                                             decreasing = T)[1:length(g)]])
write.csv(temp,file = paste0("H:/PCD_database/result4/",cancer_name,"/","node_name.csv"))

#######################要输出node 基因
gene_set_enrich_kegg_temp <-  gene_set_enrich_kegg[gene_set_enrich_kegg$Description%in%
                                                     temp ,]
kegg_gene <- gene_set_enrich_kegg_temp[,10] %>%strsplit(split = "/")
kegg_gene_1 <- list()
for(i in 1:length(kegg_gene)){
  temp <-  bitr(kegg_gene[[i]],fromType ="ENTREZID",toType = "SYMBOL",
                OrgDb = "org.Hs.eg.db" )
  kegg_gene_1[[i]] <- temp
}

names(kegg_gene_1)  <- gene_set_enrich_kegg_temp$Description
kegg_gene_1 <- do.call(rbind,kegg_gene_1)
write.csv(kegg_gene_1,paste0("H:/PCD_database/result4/",cancer_name,"/the_node_gene.csv"))  



