#########结果1结果统计和处理--------
rm(list = ls())
library(magick)
library(ggplot2)
library(tidyverse)
######### 首先对PDF图转为png
######### 图片地址
cancer_names <- list.files("H:/PCD_database/result1/")
#########癌症地址
cancer_filer <- paste0("H:/PCD_database/result1/",cancer_names)
#########癌症中图片的地址
pdf_cancer_files <- list()
for(i in 1:30){
  print(cancer_names[i])
  temp <- list.files(cancer_filer[i])
  temp <- temp[grep(".pdf",x = temp)]
  pdf_cancer_files[[i]] <- temp
}
names(pdf_cancer_files) <- cancer_names
pdf_cancer_files[[1]][1]
#####读入PDF，输出PNG 

for(i in 1:30){
  print(cancer_names[i])
  dir.create(paste0("H:/PCD_database/result1_database/",cancer_names[i]))
  for(j in 1:length(pdf_cancer_files[[i]])){
    print(pdf_cancer_files[[i]][j])
    temp_file <- paste0("H:/PCD_database/result1/",cancer_names[i],"/",pdf_cancer_files[[i]][j])
    p = image_read_pdf(path = temp_file)
    out_png <- gsub(".pdf",".png",pdf_cancer_files[[i]][j])
    temp_file_out <- paste0("H:/PCD_database/result1_database/",cancer_names[i],"/",out_png)
    image_write(p,path =temp_file_out)
  }
}

################### hallmarker的表格整理
cancer_names <- list.files("H:/PCD_database/result1/")
#########癌症地址
cancer_filer <- paste0("H:/PCD_database/result1/",cancer_names)
cancer_hall_files <- paste0(cancer_filer,"/","GSEA_HM.csv")

gesa_ha_list <- list()
for(i in 1:30){
  temp <- read.csv(cancer_hall_files[i],header = T,row.names = 1)
  gesa_ha_list[[i]] <- temp
}
gesa_ha_list <- do.call(rbind,gesa_ha_list)



gesa_ha_list_1 <- data.frame(cancer = gesa_ha_list$cancer,
                             Death  = gesa_ha_list$Cell,
                             Description = gesa_ha_list$Description,
                             enrichmentScore = gesa_ha_list$enrichmentScore,
                             NES = gesa_ha_list$NES,
                             p.adjust = gesa_ha_list$p.adjust)

write.csv(gesa_ha_list_1,
          file = "H:/PCD_database/result1_database/the_browse_cancer_pcd.csv")

#################癌症的死亡状态统计;
rm(list=ls())
library(Seurat)
cancer_names <- list.files("H:/PCD_database/result1/")
cancer_data_file <- paste0("H:/PCD_database/result1/",cancer_names,"/","PCD_ANN.RData")
cancer_death <- data.frame(cancer = 0,death_state = 0)
for(i in 1:30){
  print(i)
  load(cancer_data_file[i])
  temp <- data.frame(cancer = cancer_names[i],death_state = unique(Idents(new_data)))
  cancer_death <- rbind(cancer_death,temp)
  rm(new_data)
  gc()
}
cancer_death <- cancer_death[-1,]
write.csv(cancer_death,file = "H:/PCD_database/result1_database/cancer_death_death_state.csv")
cancer_death$umap_plot <- paste0(cancer_death$death_state,"_umap_plot.png")
cancer_death$tsne_plot <- paste0(cancer_death$death_state,"_tsne_plot.png")
cancer_death$all_tsne_plot <- "tsne_plot.png"
cancer_death$all_umap_plot <- "umap_plot.png"
cancer_death$kegg_plot <- NA

kegg <- c("Intrinsic_apoptosis",
"Extrinsic_apoptosis",
"Extrinsic_Intrinsic_apoptosis",
"Autophagy","Ferroptosis",
"Pyroptosis","Necroptosis")

for(i in 1:dim(cancer_death)[1]){
  if(cancer_death$death_state[i] =="Intrinsic_apoptosis" |
     cancer_death$death_state[i] =="Extrinsic_apoptosis" |
     cancer_death$death_state[i] =="Extrinsic_Intrinsic_apoptosis" ){
    cancer_death$kegg_plot[i] <- paste0("hsa04210.",cancer_death$cancer[i],"_",
                                     cancer_death$death_state[i],"_KEGG.png")
  }
  if(cancer_death$death_state[i] == "Autophagy" ){
    cancer_death$kegg_plot[i] <- paste0("hsa04140.",cancer_death$cancer[i],"_",
                                     cancer_death$death_state[i],"_KEGG.png")
  }
  if(cancer_death$death_state[i] == "Ferroptosis" ){
    cancer_death$kegg_plot[i] <- paste0("hsa04216.",cancer_death$cancer[i],"_",
                                        cancer_death$death_state[i],"_KEGG.png")
  }
  if(cancer_death$death_state[i] == "Pyroptosis" ){
    cancer_death$kegg_plot[i] <- paste0("hsa04621.",cancer_death$cancer[i],"_",
                                        cancer_death$death_state[i],"_KEGG.png")
  }
  if(cancer_death$death_state[i] == "Necroptosis" ){
    cancer_death$kegg_plot[i] <- paste0("hsa04217.",cancer_death$cancer[i],"_",
                                        cancer_death$death_state[i],"_KEGG.png")
  }
}
write.csv(cancer_death,file = "H:/PCD_database/result1_database/cancer_pcd_plot.csv")

################### GSEA 的一个统计图：

##########画热图
data_GSEA <- read.csv( "H:/PCD_database/result1_database/the_browse_cancer_pcd.csv")
data_GSEA <- data_GSEA[,-1]
data_GSEA$cancer_death <- paste0(data_GSEA$cancer,"_",data_GSEA$Death)
colnames(data_GSEA)
#热图绘制
matrix_gsea <- matrix(0,nrow = length(unique(data_GSEA$Description)),
                      ncol = length(unique(data_GSEA$cancer_death)))

rownames(matrix_gsea) <- unique(data_GSEA$Description)
colnames(matrix_gsea) <- unique(data_GSEA$cancer_death)

for(i in 1:dim(data_GSEA)[1]){
  matrix_gsea[data_GSEA[i,3],data_GSEA[i,7]] <- data_GSEA$NES[i]
}


library(stringr)

labels_col <- matrix(unlist(str_split(colnames(matrix_gsea),pattern = "_",
                                      n = 2)),ncol = 2,byrow = T)


annotation_col = data.frame(
  group =labels_col[,2],
  row.names = colnames(matrix_gsea))
unique(annotation_col$group)

annotation_col$group[which(annotation_col$group =="TNBC_Immunogenic_cell_death" )] <- "Immunogenic_cell_death" 
annotation_col$group[which(annotation_col$group =="TNBC_Pyroptosis" )] <- "Pyroptosis" 


ann_colors = list(
  group = c(Alkaliptosis="#41AB5D",Anoikis = "#A35D34",Autophagy="#ABCADE",
            Cuproptosis="#8DD3C7",
            Entotie_cell_death="#B3DE69",Extrinsic_apoptosis="#3B77AC",
            Ferroptosis="#BEBADA",
            Immunogenic_cell_death="#C6DBEF",Intrinsic_apoptosis="#80B1D3",
            Lysosome_dependent_cell_death="#c06f98",
            Necroptosis="#FB8072",Netotic_cell_death="#EF9BB5",
            Oxeiptosis="#F3B161",
            Parthanatos="#B7DA90",Pyroptosis = "#EB9B98",
            Extrinsic_Intrinsic_apoptosis="#CF352B",
            None="#8A959B"))

table(labels_col[,1])

c(2,6,)
p1 <- pheatmap::pheatmap(matrix_gsea,
                   cluster_cols = F,scale = "none",
                   color = colorRampPalette(c("navy","white","firebrick3"))(100),
                   border_color = "black",border=T,
                   gaps_col = c(2,6,11,15,23,25,33,42,46,50,55,62,
                                69,75,79,81,87,90,92,100,107,110,115,
                                121,124,127,129,135,140),
                   cutree_rows = 4,labels_col = labels_col[,2],
                   annotation_col=annotation_col,annotation_colors =ann_colors )
print(p1)
####################画癌症的基因
#####画marker 基因表达图
rm(list=ls())


##########载入R语言包
library(Seurat)
library(ggplot2)
library(tidyverse)
library(clusterProfiler)
##########首先需要读入PCD的程序性死亡基因
rm(list=ls())

load("D:/R-daoma/pcd_databsee/PCD_gene_set.RData")
ls()
Extrinsic_Intrinsic_apoptosis <- unique(c(Extrinsic_apoptosis,Intrinsic_apoptosis))

marker_gene_list <- list(Alkaliptosis,Anoikis,Autophagy,
                         Cuproptosis,Entotie_cell_death,Extrinsic_apoptosis,
                         Extrinsic_Intrinsic_apoptosis,Ferroptosis,Immunogenic_cell_death,
                         Intrinsic_apoptosis,Lysosome_dependent_cell_death,Necroptosis,
                         Netotic_cell_death,Oxeiptosis,Parthanatos,Pyroptosis)
aaa <- ls()
aaa <- aaa[-12]
names(marker_gene_list) <-aaa

#########读取每一个癌症的数据
file.name <- list.files("H:/PCD_database/result1/")
cancer_name  <- file.name
file.name <- paste0("H:/PCD_database/result1/",file.name,"/","PCD_ANN.RData")

data_ggplot2_point_expression <- data.frame(gene=0,data=0,group=0,cancer=0)


for(z in 1:30){
  #########读取一个癌症的数据
  print(z)
  load(file = file.name[z])
  ########提取基因表达数据
  ########按分组求平均值
  new_data <- AddMetaData(new_data,
                          metadata = Idents(new_data),
                          col.name = "Idents")
  print(cancer_name[z])
  the_death_c <- unique(new_data$Idents)
  if("None" %in%the_death_c ){
    the_death_c <- the_death_c[-which(the_death_c=="None")]}
  the_mean_data <-data.frame(gene=0,data=0,group=0)
  for(i in the_death_c){
    temp <- new_data[,colnames(new_data)[which(new_data$Idents==i)]]
    temp <- as.data.frame(temp@assays$RNA@data)
    temp1 <- temp[marker_gene_list[[i]],]
    temp1 <- na.omit(temp1)
    temp2 <- rowMeans(temp1)
    temp3 <- names(temp2[order(temp2,decreasing = T)[1:40]])
    temp4 <- data.frame(gene=temp3,data=temp2[order(temp2,decreasing = T)[1:40]],
                        group=i)
    the_mean_data <- rbind(the_mean_data,temp4)
  }
  the_mean_data <- the_mean_data[-1,]
  the_mean_data$cancer <-  cancer_name[z]
  data_ggplot2_point_expression <- rbind(data_ggplot2_point_expression,the_mean_data)
}

data_ggplot2_point_expression <- data_ggplot2_point_expression[-1,]


x <- split(data_ggplot2_point_expression,data_ggplot2_point_expression$group)
geneList <- list()
for(i in 1:length(x)){
  x_1 <-  split(x[[i]],x[[i]]$cancer)
  x_2 <- list()
  for(z in 1:length(x_1)){
    x_2[[z]] <- x_1[[z]][,1]}
  geneList[[i]] <- Reduce(intersect, x_2)
}

names(geneList) <- names(x)
which(data_ggplot2_point_expression$gene %in% unlist(geneList)  == "TRUE")


data_ggplot2_point_expression_final <- data_ggplot2_point_expression[which(data_ggplot2_point_expression$gene %in% unlist(geneList)  == "TRUE"),]

data_ggplot2_point_expression_final$gene <- factor(data_ggplot2_point_expression_final$gene,
                                                   levels = unique(unlist(geneList) ))
data_ggplot2_point_expression_final$cancer <- factor(data_ggplot2_point_expression_final$cancer,
                                                     levels = cancer_name)

data_ggplot2_point_expression_final$cancer_death <- paste0(data_ggplot2_point_expression_final$cancer,"_",
                                                           data_ggplot2_point_expression_final$group)
data_ggplot2_point_expression_final_2 <- data_ggplot2_point_expression_final[order(data_ggplot2_point_expression_final$group),]
data_ggplot2_point_expression_final$cancer_death <- factor(data_ggplot2_point_expression_final$cancer_death,
                                                           levels = unique(data_ggplot2_point_expression_final_2$cancer_death))



##############画图：
p1 <- ggplot(data_ggplot2_point_expression_final,
       aes(x=cancer_death,y = gene,size=data,
           color=data,stroke=0.2))+
  geom_point()+theme(axis.text.x=element_text(angle=90,vjust = -0.01),
                     panel.grid = element_blank(),
                     panel.background =element_blank(),
                     axis.line=element_line(colour="black", size =0.8))+
  scale_size(range = c(1,5))+
  scale_color_gradientn(colors = colorRampPalette(c("#6959CD", "#FFE7BA", 
                                                    "#B22222"))(5))+
  geom_hline(yintercept = c(2.5,12.5,17.5,30.5,36.5))+
  geom_hline(yintercept = c(8.5,9.5),linetype="dotted")+
  geom_vline(xintercept = c(18.5,26.5,49.5,57.5,67.5,95.5))

save(p1,data_ggplot2_point_expression_final,geneList,
     file = "结果一表达.RData")


#########################
######## 画每一个程序性死亡基因的气泡图
####数据读入#### 
gc()

load("D:/R-daoma/pcd_databsee/PCD_gene_set.RData")
Extrinsic_Intrinsic_apoptosis <- unique(c(Extrinsic_apoptosis,Intrinsic_apoptosis))

marker_gene_list <- list(Alkaliptosis,Anoikis,Autophagy,
                         Cuproptosis,Entotie_cell_death,Extrinsic_apoptosis,
                         Extrinsic_Intrinsic_apoptosis,Ferroptosis,Immunogenic_cell_death,
                         Intrinsic_apoptosis,Lysosome_dependent_cell_death,Necroptosis,
                         Netotic_cell_death,Oxeiptosis,Parthanatos,Pyroptosis)
aaa <- ls()
aaa <- aaa[-12]
names(marker_gene_list) <-aaa
######################计算pecrentage
i <- 16
death_state <- names(marker_gene_list)[i]
cancer_naems <- list.files(path = "H:/PCD_database/result1/")
file_data <- paste0("H:/PCD_database/result1/",cancer_naems,"/PCD_ANN.RData")
the_cancer_gene_expresion <-list()

for(x in 1:30 ){
  genes  <- marker_gene_list[[death_state]]
  print(x)
  load(file = file_data[x] )
  DefaultAssay(new_data) <-"RNA"
  new_data <-GetAssayData(object = new_data, slot = "data")[rownames(new_data)%in%genes,]
  the222 <-c()
  temp <- c()
  if(class(new_data) =="numeric"){
    load(file = file_data[x] )
    dddd <- rownames(new_data)[rownames(new_data)%in%genes]
    new_data <-GetAssayData(object = new_data, slot = "data")[rownames(new_data)%in%genes,]
    new_data <- t(as.data.frame(new_data))
    rownames(new_data) <- dddd
  }else if(dim(new_data)[1]==0){
    temp <- data.frame(Mean_expresion = 0,Percentage=0,Cancer= cancer_naems[x])
    the_cancer_gene_expresion[[x]]<- temp
    next()
  }
  
  for(j in 1:dim(new_data)[1]){
    temp <- c(temp,mean(new_data[j,]))
    the222 <- c(the222,as.numeric(mean(as.numeric(as.character(new_data[j,]))/as.numeric(as.character(colSums(new_data))),na.rm = T)))}
  temp <- data.frame(Mean_expresion = temp,Percentage=the222,Cancer= cancer_naems[x])
  rownames(temp) <- rownames(new_data)
  the_cancer_gene_expresion[[x]]<- temp
  gc()
}
save(the_cancer_gene_expresion,file = paste0(names(marker_gene_list)[i],"_pecrentage.RData"))













######################################
################细胞数量————————————
rm(list=ls())
library(Seurat)
library(plyr)
cancer_naems <- list.files(path = "H:/PCD_database/result1/")

file_data <- paste0("H:/PCD_database/result1/",cancer_naems,"/PCD_ANN.RData")


the_cell_number <- list()
for(i in 1:30){
  print(cancer_naems[i])
  load(file_data[i])
  temp <- table(Idents(new_data))/dim(new_data)[2]*100
  temp <- as.data.frame(t(as.data.frame(temp)))
  
  colnames(temp) <- temp[1,]
  temp <- temp[-1,]
  the_cell_number[[i]] <- temp
}

the_cell_number_1 <- list()
for(i in 1:30){
  print(cancer_naems[i])
  load(file_data[i])
  temp <- table(Idents(new_data))
  temp <- as.data.frame(t(as.data.frame(temp)))
  the_cell_number_1[[i]] <- temp
}


save(the_cell_number_1,file = "the_cell_number.RData")
###########图一两图合并##########
###点图###
?rbind.fill
library(plyr)
library(stringr)
rm(list=ls())

data_files<- list.files("H:/PCD_database/result1/")

data_names <- list.files("H:/PCD_database/result1/")
###############
load("the_cell_number_percent.RData")
aaa <- do.call(rbind.fill,the_cell_number)
aaa[is.na(aaa)] <- 0
rownames(aaa) <- data_names
aaa$cancer <- rownames(aaa)
library(ggplot2)
aaa$cancer <- rownames(aaa)
colnames.aaa <- colnames(aaa)[1:16]

death_cell_percent_plot_list <- list()
for(i in 1:length(colnames.aaa)){
  print(colnames.aaa[i])
  bbb <- aaa[,c(i,17)]
  bbb[,1] <- as.numeric(bbb[,1])
  bbb <- bbb[order(bbb[,1]),]
  bbb$cancer <- factor(bbb$cancer,levels=bbb$cancer)
  colnames(bbb)[1] <- "Death_state"
  
P2 <- ggplot(bbb,mapping=aes(x=cancer ,y=Death_state))+geom_bar(stat = 'identity',fill="#FABF74",width=0.7)+
    theme(axis.text.x=element_text(angle=90, hjust=1))+
    geom_text(aes(label=round(Death_state,2)), vjust=0.5,hjust=1 ,color="black", size=3)+coord_flip()+
    scale_fill_manual(values = "#FABF74")+theme_minimal()+
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())+
    theme(axis.text=element_text(size=8))+ylab(label = paste0(colnames.aaa[i],"_cells_percent"))+
    theme(axis.line = element_line(size=1, colour = "black"))+
  scale_y_continuous(expand=expansion(add = c(0.1, 0.1)),
                     limits = c(0, 100),
                     breaks = c(0,25, 50,75,100),
                     label = c("0","25", "50","75","100"))+
  theme(axis.title.x=element_text(size=12,color = "black"),
        axis.title.y= element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=10,color = "black"))+
  theme(legend.text = element_text(size =10,color = "black"),
        legend.title =element_text(size =10,color = "black") )

death_cell_percent_plot_list[[i]] <- P2

}
names(death_cell_percent_plot_list) <- colnames.aaa
save(death_cell_percent_plot_list,file = "death_cell_percent_plot_list.Rdata")
######m 1:16
load( "death_cell_percent_plot_list.Rdata")
x <- 4

colnames.aaa[x]

paste0(colnames.aaa[x],"_pecrentage.RData")
load(file =paste0(colnames.aaa[x],"_pecrentage.RData") )
for(i in 1:length(the_cancer_gene_expresion)){
  the_cancer_gene_expresion[[i]]$gene <- rownames( the_cancer_gene_expresion[[i]])
}


new_gene_expression_1 <- do.call(rbind.fill,the_cancer_gene_expresion)



#colnames(new_gene_expression_1) <- names(the_cancer_gene_expresion)
new_gene_expression_1 <- as.matrix(new_gene_expression_1)
new_gene_expression_1[which(is.na(new_gene_expression_1))] <- as.numeric(0)

new_gene_expression_1 <- as.data.frame(new_gene_expression_1)
new_gene_expression_1 <- new_gene_expression_1[grep("CASP",new_gene_expression_1$gene),]
new_gene_expression_1[,1] <- as.numeric(new_gene_expression_1[,1])
new_gene_expression_1[,2] <- as.numeric(new_gene_expression_1[,2])

library (ggplot2)
library (reshape2)#数据转换
require(scales)#数据缩放
library(ggtree)#聚类
library(aplot)#拼图
library(dplyr)


new_gene_expression_2 <- new_gene_expression_1
colnames(new_gene_expression_2)
new_gene_expression_3 <- data.frame(Mean_expresion=0,
                                    Percentage=0,Cancer=0,
                                    gene=0)


gene_u <- unique(new_gene_expression_2$gene)
for(i in gene_u){
  temp <- new_gene_expression_2[which(new_gene_expression_2$gene ==i),]
  temp$Mean_expresion <- scale(temp$Mean_expresion)
  new_gene_expression_3 <- rbind(new_gene_expression_3,temp)
}


# cancer_names <- unique(new_gene_expression_2$Cancer)
# for(i in cancer_names){
#   temp <- new_gene_expression_2[which(new_gene_expression_2$Cancer ==i),]
#   temp$Mean_expresion <- scale(temp$Mean_expresion)
#   new_gene_expression_3 <- rbind(new_gene_expression_3,temp)
# }

new_gene_expression_3 <- new_gene_expression_2

new_gene_expression_3 <- new_gene_expression_3[new_gene_expression_3$gene%in%c("CASP8","CASP6","CASP3"),]
new_gene_expression_3$gene  <- factor( new_gene_expression_3$gene,
                                       levels = names(table(new_gene_expression_3$gene)[order(table(new_gene_expression_3$gene))]))
new_gene_expression_3$Cancer <- factor(new_gene_expression_3$Cancer,levels = levels(death_cell_percent_plot_list[[colnames.aaa[x]]][["data"]][["cancer"]]))



p1 <- ggplot(new_gene_expression_3,aes(x=gene,y=Cancer))+
  geom_point(aes(size=Percentage,color=Mean_expresion))+
  geom_point(aes(size=Percentage,color=Mean_expresion),shape=1,stroke=0.1)+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black",size=1),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black",angle = 90),
        axis.text.y = element_text(size = 11,colour = "black"))+
  scale_color_gradientn(colors =colorRampPalette(c("#277AB4","#AAD0E3","white",
                                                   "#FABF74","#f2481b"))(100))+
  theme(axis.title.y = element_blank())+
  scale_y_discrete(position = "right")+
  guides(color=guide_colourbar(title = "Scale(Mean_expresion)"))+theme(legend.position="bottom") 

library(patchwork)
p1+death_cell_percent_plot_list[[colnames.aaa[x]]]+
  plot_layout(widths = c(12,2),heights = 1,ncol = 2)


###############
rm(list= ls())
load("H:/PCD_database/结果一泛癌保存数据/the_cell_number.RData")
death_naems <- c("Pyroptosis","Ferroptosis","Autophagy","Necroptosis",
                 "Cuproptosis","Parthanatos","Entotie_cell_death",
                 "Netotic_cell_death","Lysosome_dependent_cell_death",
                 "Alkaliptosis","Immunogenic_cell_death","Intrinsic_apoptosis",
                 "Extrinsic_apoptosis","Anoikis","Extrinsic_Intrinsic_apoptosis")
###############
load("H:/PCD_database/结果一泛癌保存数据/the_cell_number.RData")



num_list<- list()

for(j in 1:15){
  load(paste0("H:/PCD_database/结果一泛癌保存数据/",
                           death_naems[j],"_pecrentage.RData"))
  num <- c()
  for(i in 1:30 ){
    if(death_naems[j]%in%the_cell_number_1[[i]][1,]){
      num <- c(num,i)}}
  num_list[[j]] <- num
  }

names(num_list) <- death_naems

the_gene_num <- list()
for(i in 1:15){
  load(paste0("H:/PCD_database/结果一泛癌保存数据/",
              death_naems[i],"_pecrentage.RData"))
  for(j in 1:30){
    the_cancer_gene_expresion[[j]]$gene <- rownames( the_cancer_gene_expresion[[j]])
  }
  a <- Reduce(rbind,the_cancer_gene_expresion[num_list[[i]]])
  a <- a[which(a$Mean_expresion>0.01),]
  aa <-as.data.frame(table(a$gene))
  aa$death <- death_naems[i]
  aa <- unique(aa)
  the_gene_num[[i]] <- aa
}


aa <- do.call(rbind,the_gene_num)
write.csv(aa,"H:/PCD_database/结果一泛癌图片/death_gene_num.csv")
###############################画堆叠柱状图
rm(list = ls())
#####
cancer_naems <- list.files("H:/PCD_database/result1/")
load("H:/PCD_database/结果一泛癌保存数据/the_cell_number.RData")
the_cell_number_freq <- list()
for(i in 1:30 ){
  temp <- t(the_cell_number_1[[i]])
  temp <- as.data.frame(temp)
  temp$cancer <- cancer_naems[i]
  temp$precent<- as.numeric(temp[,2])/sum(as.numeric(temp[,2]))
  the_cell_number_freq[[i]] <- temp
}
the_cell_number_freq_matrix <- do.call(rbind,the_cell_number_freq)

library(ggplot2)

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


###########画图
table(the_cell_number_freq_matrix$cancer)
the_cell_number_freq_matrix1 <- the_cell_number_freq_matrix[the_cell_number_freq_matrix$cancer %in%c("COAD",
                                                                                                     "CRC","ESCC",
                                                                                                     "STAD","CHOL"),]

ggplot(the_cell_number_freq_matrix1, aes( x = cancer ,y=100 * precent ,
                                         fill = Var1))+
  #geom_col和geom_bar这两条命令都可以绘制堆叠柱形图
  geom_col(position = 'stack', width = 0.8)+
  #geom_bar(position = "stack", stat = "identity", width = 0.6) 
  theme_bw()+ scale_fill_manual(values=color_PCD)+
  scale_y_continuous(expand = c(0,0))+# 调整y轴属性，使柱子与X轴坐标接触+
  theme(
    text=element_text(size=12),
    plot.title = element_text(hjust = 0.5,vjust = 0.5), 
    axis.text.y=element_text(size=12,color = "black"),
    axis.text.x=element_text(size=12,  color = "black",
                             angle = 90, hjust = 0.5,
                             vjust = 0.5),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12))+
  theme_classic()

table(the_cell_number_freq_matrix1$Var1)

################################读入数据
rm(list=ls())
cancer_names <- list.files("H:/PCD_database/result1/")
file_data <- paste0("H:/PCD_database/result1/",cancer_names,"/ALL.RData")

chaojihe_enrich_AAAAAAAAA_list <- list()
ssGSEA_enrich_AAAAAAAAA_list <- list()

for(jjjjjjjjj in 1:30){
  print(jjjjjjjjj)
  load(file_data[jjjjjjjjj])
  rm(new_data)
  chaojihe_enrich_p <- as.data.frame(chaojihe_enrich_p)
  chaojihe_enrich_p$seurat_cluster <- rownames(chaojihe_enrich_p)
  chaojihe_enrich_p$cancer <- cancer_names[jjjjjjjjj]
  chaojihe_enrich_p$zhushi <- the_ident_1
  aaaaaa$seurat_cluster <- rownames(aaaaaa)
  aaaaaa$cancer <- cancer_names[jjjjjjjjj]
  aaaaaa$zhushi <- the_ident_1
  chaojihe_enrich_AAAAAAAAA_list[[jjjjjjjjj]] <- chaojihe_enrich_p
  ssGSEA_enrich_AAAAAAAAA_list[[jjjjjjjjj]] <- aaaaaa
}
save(chaojihe_enrich_AAAAAAAAA_list,ssGSEA_enrich_AAAAAAAAA_list,file = "数据库结果1.RData")

load("数据库结果1.RData")
library(plyr)
chaojihe_enrich_AAAAAAAAA_list_111111 <- do.call(rbind.fill,chaojihe_enrich_AAAAAAAAA_list)
chaojihe_enrich_AAAAAAAAA_list_111111 <- chaojihe_enrich_AAAAAAAAA_list_111111[,-16]
ssGSEA_enrich_AAAAAAAAA_list_111111 <- do.call(rbind.fill,ssGSEA_enrich_AAAAAAAAA_list)
ssGSEA_enrich_AAAAAAAAA_list_111111 <- ssGSEA_enrich_AAAAAAAAA_list_111111[,-16]
library(tidyverse)
chaojihe_enrich_mean <-  chaojihe_enrich_AAAAAAAAA_list_111111 %>% group_by(cancer,zhushi) %>% summarise_all(mean)
ssGSEA_enrich_mean <- ssGSEA_enrich_AAAAAAAAA_list_111111 %>% group_by(cancer,zhushi) %>% summarise_all(mean)
library(ggplot2)
##########画一种
load("death_cell_percent_plot_list.Rdata")




########Pyroptosis 
for(kkkkkkkk in 1:14){
death_names <- names(death_cell_percent_plot_list)
death_names <- death_names[2:16]
kkkkkkk<- kkkkkkkk
death_names[kkkkkkk]
temp1 <- chaojihe_enrich_mean%>%filter(zhushi ==death_names[kkkkkkk])
temp1 <- temp1[,c("cancer",death_names[kkkkkkk])]
temp2 <- data_frame(cancer = unique(chaojihe_enrich_mean$cancer)[!unique(chaojihe_enrich_mean$cancer)%in%
                                                                   temp1$cancer],p = 1)
colnames(temp2)[2] <- death_names[kkkkkkk]
temp3 <- rbind(temp1,temp2)
temp3$x= 1
temp3$cancer <- factor(temp3$cancer,levels = levels(death_cell_percent_plot_list[[death_names[kkkkkkk]]][["data"]][["cancer"]]))
colnames(temp3)[2] <- "Immunogenic_cell_death"
temp3$Immunogenic_cell_death[which(temp3$Immunogenic_cell_death==0)] <- 0.0000000001
p1 <- ggplot(temp3,aes(x=x,y=cancer))+
  geom_point(aes(size=-log10(Immunogenic_cell_death),
                 color=-log10(Immunogenic_cell_death)))+ 
  xlab(death_names[kkkkkkk])+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black",size=1),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x.bottom = element_text(size = 11,colour = "black",angle = 90),
        axis.text.y = element_text(size = 11,colour = "black"))+
  scale_color_gradientn(colors =colorRampPalette(c("grey","#AAD0E3","#277AB4",
                                                   "#FABF74","#f2481b"))(100))+
  theme(axis.title.y = element_blank())+
  scale_y_discrete(position = "right")+
  guides(color=guide_colourbar(title = "-log10(p.adjust)"),
         size=guide_legend(title = "-log10(p.adjust)"))+
  theme(legend.position="left") 
  
temp4 <- ssGSEA_enrich_mean%>%filter(zhushi ==death_names[kkkkkkk])
temp4 <- temp4[,c("cancer",death_names[kkkkkkk])]
temp5 <- data_frame(cancer = unique(chaojihe_enrich_mean$cancer)[!unique(chaojihe_enrich_mean$cancer)%in%
                                                                   temp1$cancer],NES = 0)
colnames(temp5)[2] <- death_names[kkkkkkk]
temp6 <- rbind(temp4,temp5)
temp6$x= as.numeric(1)
colnames(temp6)[2] <- "Immunogenic_cell_death"
temp6$cancer <- factor(temp6$cancer,levels = levels(death_cell_percent_plot_list[[death_names[kkkkkkk]]][["data"]][["cancer"]]))
#temp6$Immunogenic_cell_death[which(temp6$Immunogenic_cell_death==0)] <- 0.0000000001
p2 <- ggplot(temp6,aes(x=x,y=cancer))+
  geom_point(aes(size=abs(Immunogenic_cell_death),
                 color=abs(Immunogenic_cell_death)))+ 
  xlab(death_names[kkkkkkk])+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black",size=1),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x.bottom = element_text(size = 11,colour = "black",angle = 90),
        axis.text.y = element_text(size = 11,colour = "black"))+
  scale_color_gradientn(colors =colorRampPalette(c("grey","#AAD0E3",
                                                   "#FABF74","#fc00ff"))(100))+
  theme(axis.title.y = element_blank())+
  scale_y_discrete(position = "left")+
  guides(color=guide_colourbar(title = "abs(ssGSEA_score)"),
         size=guide_legend(title = "abs(ssGSEA_score)"))+
  theme(legend.position="right")


library(patchwork)
p5 <- p2+p1+death_cell_percent_plot_list[[death_names[kkkkkkk]]]+
  plot_layout(widths = c(2,2,8),heights = 1,ncol = 3)

pdf(file = paste0(death_names[kkkkkkk],".pdf"),
    width = 11.69,height = 8.27)
print(p5)
dev.off()
}

#########################   Extrinsic_Intrinsic_apoptosis_pecrentage  这个要两个一起

kkkkkkk<- 15
death_names[kkkkkkk]
temp1 <- chaojihe_enrich_mean%>%filter(zhushi %in% c("Extrinsic_apoptosis","Intrinsic_apoptosis"))
temp1 <- temp1[,c("cancer", c("Extrinsic_apoptosis","Intrinsic_apoptosis"))]

temp1$mean <- rowMeans(temp1[,2:3])
temp1 <- temp1[,c(1,4)]
temp2 <- data_frame(cancer = unique(chaojihe_enrich_mean$cancer)[!unique(chaojihe_enrich_mean$cancer)%in%
                                                                   temp1$cancer],p = 1)
colnames(temp2)[2] <- "mean"
temp3 <- rbind(temp1,temp2)
temp3$x= 1
temp3$cancer <- factor(temp3$cancer,levels = levels(death_cell_percent_plot_list[[death_names[kkkkkkk]]][["data"]][["cancer"]]))
colnames(temp3)[2] <- "Immunogenic_cell_death"
temp3$Immunogenic_cell_death[which(temp3$Immunogenic_cell_death==0)] <- 0.0000000001
temp3_1 <-  temp3 %>% group_by(cancer) %>% summarise_all(mean)


p1 <- ggplot(temp3_1,aes(x=x,y=cancer))+
  geom_point(aes(size=-log10(Immunogenic_cell_death),
                 color=-log10(Immunogenic_cell_death)))+ 
  xlab(death_names[kkkkkkk])+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black",size=1),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x.bottom = element_text(size = 11,colour = "black",angle = 90),
        axis.text.y = element_text(size = 11,colour = "black"))+
  scale_color_gradientn(colors =colorRampPalette(c("grey","#AAD0E3","#277AB4",
                                                   "#FABF74","#f2481b"))(100))+
  theme(axis.title.y = element_blank())+
  scale_y_discrete(position = "right")+
  guides(color=guide_colourbar(title = "-log10(p.adjust)"),
         size=guide_legend(title = "-log10(p.adjust)"))+
  theme(legend.position="left") 

temp4 <- ssGSEA_enrich_mean%>%filter(zhushi %in% c("Extrinsic_apoptosis","Intrinsic_apoptosis"))
temp4 <- temp4[,c("cancer",c("Extrinsic_apoptosis","Intrinsic_apoptosis"))]
temp4$mean <- rowMeans(temp4[,c(2,3)])
temp4 <- temp4[,c(1,4)]
temp5 <- data_frame(cancer = unique(chaojihe_enrich_mean$cancer)[!unique(chaojihe_enrich_mean$cancer)%in%
                                                                   temp1$cancer],NES = 0)
colnames(temp5)[2] <-"mean"
temp6 <- rbind(temp4,temp5)
temp6$x= as.numeric(1)
colnames(temp6)[2] <- "Immunogenic_cell_death"
temp6$cancer <- factor(temp6$cancer,levels = levels(death_cell_percent_plot_list[[death_names[kkkkkkk]]][["data"]][["cancer"]]))
temp6_1 <-  temp6 %>% group_by(cancer) %>% summarise_all(mean)

#temp6$Immunogenic_cell_death[which(temp6$Immunogenic_cell_death==0)] <- 0.0000000001
p2 <- ggplot(temp6,aes(x=x,y=cancer))+
  geom_point(aes(size=abs(Immunogenic_cell_death),
                 color=abs(Immunogenic_cell_death)))+ 
  xlab(death_names[kkkkkkk])+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black",size=1),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x.bottom = element_text(size = 11,colour = "black",angle = 90),
        axis.text.y = element_text(size = 11,colour = "black"))+
  scale_color_gradientn(colors =colorRampPalette(c("grey","#AAD0E3",
                                                   "#FABF74","#fc00ff"))(100))+
  theme(axis.title.y = element_blank())+
  scale_y_discrete(position = "left")+
  guides(color=guide_colourbar(title = "abs(ssGSEA_score)"),
         size=guide_legend(title = "abs(ssGSEA_score)"))+
  theme(legend.position="right")



library(patchwork)
p2+p1+death_cell_percent_plot_list[[death_names[kkkkkkk]]]+
  plot_layout(widths = c(2,2,8),heights = 1,ncol = 3)


##################画动态柱状图：
require(devtools)
chooseCRANmirror()
BiocManager::install("palmerpenguins")
install_github('yihui/recharts')
install.packages('plotly') 
install_github('Lchiffon/REmap')
install.packages('igraph') 
install.packages('networkD3')
BiocManager::install("htmlwidgets")
library(htmlwidgets)
library(plotly)

#aa<-aggregate(paste(colnc$lncRNA1,colnc$lncRNA2,sep = "_"),list(colnc$cell_type,colnc$tumor),length)

bb<-reshape2::dcast(the_cell_number_freq_matrix,
                    cancer~Var1,value.var = "Freq")
bb[,2:17] <- apply(bb[,2:17],c(1,2),as.numeric)
bb[is.na(bb)] <- 0
#cell<-intersect(unique(the_cell_number_freq_matrix$cancer),unique(aa$Group.1))
cell<-unique(the_cell_number_freq_matrix$cancer)

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
fig<-plot_ly()

aaaaaa <- bb[,c(1,3)]
aaaaaa$text <- paste("Cancer: ", aaaaaa[,1],
  "<br>Death State: ", aaaaaa[,2],
  "<br>Cell Number: ", aaaaaa[,2])

p <- plot_ly()

p <- p %>% add_trace(x = aaaaaa$cancer,y = aaaaaa[,2],
                  type = "bar", name = colnames(aaaaaa)[2],text = aaaaaa$text,
                  marker = list(color = color_PCD[colnames(aaaaaa)[2]])) %>% 
  layout( xaxis = list(title = ""),yaxis = list(title = ""))

p <- p %>% add_trace(x = aaaaaa$cancer,y = aaaaaa[,2],
                   type = "bar", name = colnames(aaaaaa)[2],text = aaaaaa$text,
                   marker = list(color = color_PCD[colnames(aaaaaa)[2]]))
?plot_ly
p <- plot_ly()

for(i in 1:16){
  aaaaaa <- bb[,c(1,i+1)]
  aaaaaa$text <- paste("Cancer: ", aaaaaa[,1],
                       "<br>Death State: ", colnames(aaaaaa)[2],
                       "<br>Cell Number: ", aaaaaa[,2])
  p <- p %>% add_trace(x = aaaaaa$cancer,y = aaaaaa[,2],
                       type = "bar", name = colnames(aaaaaa)[2],
                       text = aaaaaa$text,
                       marker = list(color = color_PCD[colnames(aaaaaa)[2]]),
                       textposition ="none") %>% 
    layout( xaxis = list(title = ""),yaxis = list(title = ""))
  
}



p <- p %>% layout(yaxis = list(title = 'Count'),xaxis = list(title = "Cancer"), barmode = 'stack')
p <- p %>% layout(legend = list(y = 1.3,orientation = 'h'))
p

#chooseCRANmirror()
saveWidget(partial_bundle(p), "主页动态统计.html")
####################画每个癌症
the_cell_number_freq_matrix 
cancer_names <- list.files("H:/PCD_database/result1/")
for(i in 1:30){
  
temp <- the_cell_number_freq_matrix[the_cell_number_freq_matrix$cancer %in% cancer_names[i] ,]


temp$text <- paste("Cancer: ", temp$cancer,
                     "<br>Death State: ",temp$Var1 ,
                     "<br>Cell Number: ", temp$Freq)
temp$Freq <- as.numeric(temp$Freq)
temp$color <- color_PCD[temp$Var1]
temp <- temp[order(temp$Var1,decreasing = T),]
temp$Var1 <- factor(temp$Var1,levels = temp$Var1)
temp$color <- factor(temp$color,levels = temp$color)
color <-  temp$color
names(color) <-temp$Var1
p <- plot_ly()
p <- p %>% add_trace(x = temp$Freq,y = temp$cancer,
                     type = "bar", name = temp$Var1,
                     text = temp$text,
                     marker = list(color = color),
                     textposition ="none") %>% 
  layout( xaxis = list(title = ""),yaxis = list(title = ""))

p <- p %>% layout(yaxis = list(title = 'Cancer'),
                  xaxis = list(title = "Count"), 
                  barmode = 'stack')
p <- p %>% layout(legend = list(y = 1.3,orientation = 'h'))

saveWidget(partial_bundle(p), paste0("H:/PCD_database/动态图/每个癌症/",cancer_names[i],"统计.html"))
}

rm(list=ls())
library(htmlwidgets)
library(plotly)
########画每一个死亡特征在癌症中的比例图：
########首先需要读入数据
load("death_cell_percent_plot_list.Rdata")
load("the_cell_number_percent.RData")
load("the_cell_number.RData")
cancer_naems <- as.character(death_cell_percent_plot_list[["Extrinsic_apoptosis"]][["data"]][["cancer"]])
cancer_naems <- sort(cancer_naems)
the_cell_number_freq <- list()
for(i in 1:30 ){
  temp <- t(the_cell_number_1[[i]])
  temp <- as.data.frame(temp)
  temp$cancer <- cancer_naems[i]
  temp$precent<- as.numeric(temp[,2])/sum(as.numeric(temp[,2]))
  the_cell_number_freq[[i]] <- temp
}
the_cell_number_freq_matrix <- do.call(rbind,the_cell_number_freq)
sum(as.numeric(the_cell_number_freq_matrix$Freq))
############## 

bb<-reshape2::dcast(the_cell_number_freq_matrix,
                    cancer~Var1,value.var = "precent")

bb[,2:17] <- apply(bb[,2:17],c(1,2),as.numeric)
bb[is.na(bb)] <- 0

dd <-reshape2::dcast(the_cell_number_freq_matrix,
                     cancer~Var1,value.var = "Freq")
dd[,2:17] <- apply(dd[,2:17],c(1,2),as.numeric)
dd[is.na(dd)] <- 0


#cell<-intersect(unique(the_cell_number_freq_matrix$cancer),unique(aa$Group.1))
cell<-unique(the_cell_number_freq_matrix$cancer)

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
# fig<-plot_ly()
# 
# aaaaaa <- bb[,c(1,3)]
# aaaaaa$text <- paste("Cancer: ", aaaaaa[,1],
#                      "<br>Death State: ", aaaaaa[,2],
#                      "<br>Cell Number: ", aaaaaa[,2])
# 
# p <- plot_ly()
# 
# p <- p %>% add_trace(x = aaaaaa$cancer,y = aaaaaa[,2],
#                      type = "bar", name = colnames(aaaaaa)[2],text = aaaaaa$text,
#                      marker = list(color = color_PCD[colnames(aaaaaa)[2]])) %>% 
#   layout( xaxis = list(title = ""),yaxis = list(title = ""))
# 
# p <- p %>% add_trace(x = aaaaaa$cancer,y = aaaaaa[,2],
#                      type = "bar", name = colnames(aaaaaa)[2],text = aaaaaa$text,
#                      marker = list(color = color_PCD[colnames(aaaaaa)[2]]))
# ?plot_ly


########选定几个颜色
0-500"grey"
1——1000   "#6295A2"
1001-2000 "#5A639C"
2001-5000 "#FFAF45"
>5000.   "#D24545"

for(i in 1:16){
  print(i)
  aaaaaa <- bb[,c(1,i+1)]
  aaaaaa$num <- dd[,i+1]
  aaaaaa$text <- paste("Cancer: ", aaaaaa[,1],
                       "<br>Death State: ", colnames(aaaaaa)[2],
                       "<br>Cell Number: ", aaaaaa[,3],
                       "<br>Cell Precent: ", aaaaaa[,2])
  aaaaaa <- aaaaaa[-which(aaaaaa[,2]==0),]
  temp  <- aaaaaa[order(aaaaaa[,2],decreasing = T),]
  temp$cancer <- factor(temp$cancer,levels = temp$cancer)
  temp$color <- NA
  temp$name <- NA
  for(j in 1:dim(temp)[1]){
    if(temp$num[j] <= 500){
      temp$color[j] <-"grey"
      temp$name[j] <- "(0,500]"
    } else if(temp$num[j] > 500 & temp$num[j] <= 1000){
      temp$color[j] <-"#6295A2"
      temp$name[j] <- "(500,1000]"
    }else if(temp$num[j] > 1000 & temp$num[j] <= 2000){
      temp$color[j] <-"#5A639C"
      temp$name[j] <- "(1000,2000]"
    }else if(temp$num[j] > 2000 & temp$num[j] <= 5000){
      temp$color[j] <-"#FFAF45"
      temp$name[j] <- "(2000,5000]"
    }else if(temp$num[j] > 5000){
      temp$color[j] <-"#D24545"
      temp$name[j] <- "(5000,50000]"
    }
  }
  
  color <-  c("grey",
              "#6295A2",
              "#5A639C",
              "#FFAF45",
              "#D24545")
  
  names(color) <-c("(0,500]","(500,1000]","(1000,2000]","(2000,5000]","(5000,50000]")
  
  sort(unique(temp$name))
  p <- plot_ly(x = temp$cancer,y = temp[,2]*100,
               type = "bar",colors = color[sort(unique(temp$name))],
               color = factor(temp$name,levels = c("(0,500]","(500,1000]","(1000,2000]","(2000,5000]","(5000,50000]") ),
               text = temp$text,
               textposition ="none")%>%
    layout( xaxis = list(title = "Cancer"),yaxis = list(title = "Precent(%)"),
            legend = list(y = 1.3,orientation = 'h'))
  
  p <- p %>% add_annotations(text=paste0(round(temp[,2],3)*100,"%"), showarrow=F,
                             textangle = 90,yanchor= "bottom")
  p
  saveWidget(partial_bundle(p), file = paste0("cell death precent/",colnames(temp)[2],"_precent.html"))
  # p <- plot_ly(x = temp$cancer,y = temp[,2],
  #                    type = "bar",name =  color,
  #                    text = temp$text,
  #                    marker = list(color = color),
  #                    textposition ="none")%>% 
  #   layout(title = paste0("The Precent of " ,colnames(temp)[2]," State Cells"))%>%
  #   layout( xaxis = list(title = "Cancer"),yaxis = list(title = "Precent"),
  #           legend = list(y = 1.3,orientation = 'h'))
  #   
}


death_cell_percent_file <- data.frame(death_state = colnames(bb)[2:17],file = paste0("https://secretdatabase.pages.dev/",colnames(bb)[2:17],"_precent.html") )

write.csv(death_cell_percent_file,file = "death_cell_percent_file.csv")



