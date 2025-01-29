###################结果二整理------
rm(list=ls())
###################轨迹的统计
guiji <- read.csv("H:/PCD_database/轨迹.csv")
rownames(guiji) <- guiji$cancer
rownames(guiji)[6] <- "BRCA_TNBC"
guiji$cancer[6] <- "BRCA_TNBC"
guiji$trajectory <- paste0(guiji$from,"->",guiji$to)
length(table(guiji$trajectory))
aaa <- as.data.frame(table(guiji$trajectory))
write.csv(guiji,file = "TEMP.csv")
#######pheatmap ###
the_heatmap <- matrix(ncol = 30,nrow = 22)
colnames(the_heatmap) <- guiji$cancer
rownames(the_heatmap) <- unique(guiji$trajectory)

for(i in colnames(the_heatmap)){
  the_heatmap[guiji[i,"trajectory"],i] <- as.numeric(1)
}
the_heatmap[is.na(the_heatmap)] <- as.numeric(0)

library(pheatmap)

pheatmap(t(the_heatmap),cluster_rows = T,cluster_cols = T,
         color = colorRampPalette(c("white","firebrick3"))(100)
         ,display_numbers = F,scale="none",
         show_rownames = T,
         show_colnames = T,
         border= T,
         cellwidth = 10,
         cellheight = 10)

###############共享基因#####
#####获取每个癌症轨迹的差异基因的地址
the_de_pr <- list.dirs("H:/PCD_database/result2/")
the_de_pr <- the_de_pr[-1]
the_de_pr <- paste0(the_de_pr, "/pr_graph_test_res.RData")

list_module <- list()
for (i in guiji$trajectory) {
  temp <-  rownames(guiji)[which(guiji$trajectory == i)]
  temp1 <- which(guiji$trajectory == i)
  temp3 <- c()
  list_temp <- list()
  for (j in 1:length(temp1)) {
    temp2 <- get(load(the_de_pr[temp1[j]]))
    pr_graph_test_res_1 <- temp2[- which(abs(temp2$morans_I)<0.15),]
    temp2 <- row.names(subset(pr_graph_test_res_1, q_value < 0.001))
    list_temp[[j]] <- temp2
  }
  if (length(list_temp) > 1) {
    list_module[[i]] <- Reduce(intersect, list_temp)
  } else{
    list_module[[i]] <- list_temp[[1]]
  }
}
save(list_module,file = "22种转化轨迹的癌症共享差异基因.RData")
#load("22种转化轨迹的癌症共享差异基因.RData")
######进行Go功能富集：
library(clusterProfiler)
library(org.Hs.eg.db)
gene_enrich_list <- list()
for (i in 1:22) {
  temp <- enrichGO(
    gene = list_module[[i]],
    keyType = "SYMBOL",
    OrgDb = "org.Hs.eg.db",
    ont = "BP"
  )
  temp <- as.data.frame(temp)
  temp$modules <- names(list_module)[i]
  gene_enrich_list[[i]] <- temp[1:300,]
}

the_Go_enrich_jinhua_modules <- do.call(rbind, gene_enrich_list)
the_Go_enrich_jinhua_modules <- na.omit(the_Go_enrich_jinhua_modules)
write.csv(the_Go_enrich_jinhua_modules,"进化模块的共享差异基因功能.csv")


######
library(ggplot2)
for (i in 1:22) {
  i <- 14
  p <- ggplot(data = gene_enrich_list[[i]][1:30,],
              aes(
                x = Description,
                y = Count,
                fill = "#FD8D62")) +
    geom_bar(stat = "identity", width = 0.8)  +
    scale_fill_manual(values = "#FD8D62") + theme_bw() +
    xlab("GO term") + ylab("Num of Genes") +
    labs(title = "The Most Enriched GO Terms") +
    theme(axis.text.x = element_text(
      color = "black",
      angle = 70,
      vjust = 1,
      hjust = 1),
      axis.title.y = element_text(color = "black"))+
    guides(fill=FALSE) 
  dev.new()
  pdf(file = paste0("H:/PCD_database/结果二泛癌图片/",
      gsub(
        pattern = ">" ,
        replacement = "" ,
        x = names(list_module)[i]
      ),
      "_go_enrich.pdf"
    ),
    width = 7,
    height = 7
  )
  print(p)
  dev.off()
}

load("D:/R-daoma/pcd_databsee/PCD_gene_set.RData")
marker_gene_list <- list(Pyroptosis=Pyroptosis,
                         Ferroptosis=Ferroptosis,
                         Autophagy=Autophagy,Necroptosis=Necroptosis,Cuproptosis=Cuproptosis,
                         Parthanatos=Parthanatos,
                         Entotie_cell_death=Entotie_cell_death,
                         Netotic_cell_death=Netotic_cell_death,
                         Lysosome_dependent_cell_death=Lysosome_dependent_cell_death,
                         Alkaliptosis=Alkaliptosis,Oxeiptosis=Oxeiptosis,
                         Immunogenic_cell_death=Immunogenic_cell_death,
                         Intrinsic_apoptosis=Intrinsic_apoptosis
                         ,Extrinsic_apoptosis=Extrinsic_apoptosis,
                         Anoikis=Anoikis,
                         Extrinsic_Intrinsic_apoptosis=c(unique(c(Extrinsic_apoptosis,
                                                                  Intrinsic_apoptosis))))

######基因画图

library(monocle3)
i <- 15
names(list_module)[i]
load("D:/R-daoma/pcd_databsee/PCD_gene_set.RData")

which(guiji$trajectory == names(list_module)[i])
######基因需要修改
####单个癌症
gene_temp <- intersect(list_module[[i]], marker_gene_list$Pyroptosis)
gene_temp <- gene_temp[sample(1:length(gene_temp), 12)]
####多癌症

gene_temp_1 <-intersect(list_module[[i]], marker_gene_list$Entotie_cell_death)

gene_temp_2 <-intersect(list_module[[i]], marker_gene_list$Extrinsic_apoptosis)

gene_temp <- c(gene_temp_2[sample(1:length(gene_temp_2), 10)],
               gene_temp_1[sample(1:length(gene_temp_1), 2)])

gene_temp <- c(gene_temp_1,gene_temp_2,sample(list_module[[i]],3))
gene_temp <- c("STMN1","EGR1","IRF3")
#gene_temp <- c(sample(list_module[[i]], 12))
# 
# gene_temp <- gene_temp_2[1:12]


######单个癌症
# the_module <- paste0("H:/PCD_database/result2", "/", guiji$cancer, "/CDS.RData")
# names(the_module) <-  guiji$cancer
the_module[guiji[which(guiji$trajectory == names(list_module)[i]),1]]

load(the_module[guiji[which(guiji$trajectory == names(list_module)[i]),1]])
pes_time <- pseudotime(cds)
data_temp <- SingleCellExperiment::counts(cds[gene_temp,])
data_temp  <- as.matrix(data_temp)

data_temp <- t(data_temp)
data_temp <- apply(data_temp, c(1, 2), as.numeric)
data_temp[which(data_temp == 0)] <- NA
data_temp <- as.data.frame(data_temp)
data_temp$time <- as.numeric(pes_time)
names(the_module)[which(guiji$trajectory == names(list_module)[i])]
data_temp$cancer <- names(the_module)[which(guiji$trajectory == names(list_module)[i])]
fffff <- colnames(data_temp)

p_list <- list()

for (i in 1:12) {
  print(i)
  data_temp_1 <- data_temp[, c(i, 13, 14)]
  colnames(data_temp_1)[1] <- "GENE"
  p <-
    ggplot(data = data_temp_1, aes(x = time, y = GENE, color = time)) +
    geom_point(size =2) +
    geom_smooth(aes(fill = cancer), na.rm = T,color="#11998e")+ 
    theme_classic() +
    scale_color_gradient(low = "#FFDD94", high = "#FD8D62") +
    theme(axis.text.x = element_text(size = 13,color = "black",hjust = 0.5),
      axis.text.y = element_text(size = 13,color = "black",hjust = 0.5) ,
      axis.line.x = element_line(linetype = 1,color = "black",size = 1),
      axis.line.y = element_line(linetype = 1,color = "black",size = 1)) +
    labs(title = colnames(data_temp[i]))
  p_list[[i]] <- p+scale_fill_manual(values = c("#11998e"))
  
}



library(ggpubr)
ggarrange(
  p_list[[1]],
  p_list[[2]] , 
  p_list[[3]],
  p_list[[4]],
  p_list[[5]],
  p_list[[6]],
  p_list[[7]],
  p_list[[8]],
  p_list[[9]],
  p_list[[10]],
  p_list[[11]],
  p_list[[12]],
  common.legend = TRUE
)



######多个癌症

temp_123 <- which(guiji$trajectory == names(list_module)[i])

names(the_module)[temp_123]

data_gene <- list()
for (i in temp_123) {
  load(the_module[i])
  pes_time <- pseudotime(cds)
  data_temp <- SingleCellExperiment::counts(cds[gene_temp,])
  data_temp  <- as.matrix(data_temp)
  data_temp <- t(data_temp)
  data_temp <- apply(data_temp, c(1, 2), as.numeric)
  data_temp[which(data_temp == 0)] <- NA
  data_temp <- as.data.frame(data_temp)
  data_temp$time <- as.numeric(pes_time)
  data_temp$cancer <- names(the_module)[i]
  data_gene[[names(the_module)[i]]] <- data_temp
}
gene_temp%in%rownames(cds)

data_temp <- do.call(rbind, data_gene)

fffff <- colnames(data_temp)

p_list <- list()
# ggplot(data = data_temp_1, aes(x = time, y = GENE, color = time)) +
#   geom_point(size =2) +
#   geom_smooth(aes(fill = cancer), na.rm = T)

for (i in 1:3) {
  print(i)
  data_temp_1 <- data_temp[, c(i, 4, 5)]
  colnames(data_temp_1)[1] <- "GENE"
  p <-
    ggplot(data = data_temp_1, aes(x = time, y = GENE, color = time)) +
    geom_point(size =2) +
    geom_smooth(aes(fill = cancer), 
                na.rm = T,
                color=c("#11998e"))+ 
    theme_classic() +
    scale_color_gradient(low = "#FFDD94", high = "#FD8D62") +
    theme(axis.text.x = element_text(size = 13,color = "black",hjust = 0.5),
          axis.text.y = element_text(size = 13,color = "black",hjust = 0.5) ,
          axis.line.x = element_line(linetype = 1,color = "black",size = 1),
          axis.line.y = element_line(linetype = 1,color = "black",size = 1)) +
    labs(title = colnames(data_temp[i]))+ylim(c(0,10))
  p_list[[i]] <- p+scale_fill_manual(values = c(LIHC="#283c86",OV="#45a247"
                                               ,PPB="#5C258D",PRAD="#EB3349"))
}

library(ggpubr)
ggarrange(
  p_list[[1]],
  p_list[[2]],
  p_list[[3]],
  common.legend = TRUE
)




#########################################
rm(list = ls())
guiji <- read.csv("H:/PCD_database/轨迹.csv")
rownames(guiji) <- guiji$cancer
rownames(guiji)[6] <- "BRCA_TNBC"
guiji$cancer[6] <- "BRCA_TNBC"
guiji$trajectory <- paste0(guiji$from,"->",guiji$to)

###############对泛癌结果的图片进行处理
library(magick)
trajectory <- unique(guiji$trajectory)
Go_enrich <- paste0("H:/PCD_database/结果二泛癌图片/",gsub(">","",trajectory),"_go_enrich.pdf")
gene_plot_file <- paste0("H:/PCD_database/结果二泛癌图片/",gsub(">","",trajectory),"_gene.pdf")

Go_enrich_png <- paste0("H:/PCD_database/结果二泛癌图片/",gsub(">","",trajectory),"_go_enrich.png")
gene_plot_file_png <- paste0("H:/PCD_database/结果二泛癌图片/",gsub(">","",trajectory),"_gene.png")

for(i in 1:22){
    print(trajectory[i])
    p = image_read_pdf(path = Go_enrich[i])
    image_write(p,path =Go_enrich_png[i])
    p = image_read_pdf(path = gene_plot_file[i])
    image_write(p,path =gene_plot_file_png[i])
}

plot_file <- data.frame(modules = trajectory,
                        Go_plot = paste0(gsub(">","",trajectory),"_go_enrich.png"),
                        Gene_plot = paste0(gsub(">","",trajectory),"_gene.png"))

write.csv(plot_file,"H:/PCD_database/泛癌整理/死亡转化模块/plot_file.csv")

#####################结果二的单个癌症图片
rm(list=ls())
cancer_names <- list.files("H:/PCD_database/result2/")

pseudotime_plot_pdf_files <- paste0("H:/PCD_database/result2/",cancer_names,"/pseudotime_plot.pdf")
pseudotime_plot_png_files <- paste0("H:/PCD_database/result2_database/",cancer_names,"/pseudotime_plot.png")


pcd_plot_pdf_files <- paste0("H:/PCD_database/result2/",cancer_names,"/pcd_plot.pdf")
pcd_plot_png_files <- paste0("H:/PCD_database/result2_database/",cancer_names,"/pcd_plot.png")


malignant_plot_pdf_files <- paste0("H:/PCD_database/result2/",cancer_names,"/malignant_plot.pdf")
malignant_plot_png_files <- paste0("H:/PCD_database/result2_database/",cancer_names,"/malignant_plot.png")


GO_BP_plot_1_pdf_files <- paste0("H:/PCD_database/result2/",cancer_names,"/GO_BP_plot_1.pdf")
GO_BP_plot_1_png_files <- paste0("H:/PCD_database/result2_database/",cancer_names,"/GO_BP_plot_1.png")


GO_BP_plot_2_pdf_files <- paste0("H:/PCD_database/result2/",cancer_names,"/GO_BP_plot_2.pdf")
GO_BP_plot_2_png_files <- paste0("H:/PCD_database/result2_database/",cancer_names,"/GO_BP_plot_2.png")


GO_BP_plot_3_pdf_files <- paste0("H:/PCD_database/result2/",cancer_names,"/GO_BP_plot_3.pdf")
GO_BP_plot_3_png_files <- paste0("H:/PCD_database/result2_database/",cancer_names,"/GO_BP_plot_3.png")


for(i in 1:30 ){
  print(cancer_names[i])
  dir.create(paste0("H:/PCD_database/result2_database/",cancer_names[i]))
  p = image_read_pdf(path = pseudotime_plot_pdf_files[i])
  image_write(p,path =pseudotime_plot_png_files[i])
  
  p = image_read_pdf(path = pcd_plot_pdf_files[i])
  image_write(p,path =pcd_plot_png_files[i])
  
  p = image_read_pdf(path = malignant_plot_pdf_files[i])
  image_write(p,path =malignant_plot_png_files[i])
  
  p = image_read_pdf(path = GO_BP_plot_1_pdf_files[i])
  image_write(p,path =GO_BP_plot_1_png_files[i])
  
  p = image_read_pdf(path = GO_BP_plot_2_pdf_files[i])
  image_write(p,path =GO_BP_plot_2_png_files[i])
  
  p = image_read_pdf(path = GO_BP_plot_3_pdf_files[i])
  image_write(p,path =GO_BP_plot_3_png_files[i])
}

plot_file <- data.frame(cancer=cancer_names,pseudotime_plot="pseudotime_plot.png",
                        malignant_plot = "malignant_plot.png",pcd_plot = "pcd_plot.png",
                        GO_BP_plot_1="GO_BP_plot_1.png", GO_BP_plot_2="GO_BP_plot_2.png",
                        GO_BP_plot_3="GO_BP_plot_3.png")
write.csv(plot_file ,file = "H:/PCD_database/result2_database/plot_file.csv")




##############对单个癌症进行一个整理—
rm(list=ls())
library(dplyr)
library(clusterProfiler)
cancer_names <- list.files("H:/PCD_database/result2/")
enrichmentData_list <- list()
for(i in 1:30){
  load(paste0("H:/PCD_database/result2/",cancer_names[i],"/GO_BP.RData"))
  enrichmentData <- enrich@result %>%
    as.data.table()
  enrichmentData <- enrichmentData[-which(enrichmentData$p.adjust>0.05),]
  enrichmentData$cancer <- cancer_names[i]
  enrichmentData_list[[i]] <-enrichmentData
}
enrichmentData_list_save <- do.call(rbind,enrichmentData_list)
write.csv(enrichmentData_list_save,file = "H:/PCD_database/result2_database/GO_BP.csv")











