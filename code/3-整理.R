#########结果三的数据整理
rm(list=ls())
library(magick)

####################地址文件的设置
cancer_names <- list.files("H:/PCD_database/result3/")
######################对PDF的数据进行png转化
plot_file <- data.frame(cancer = cancer_names,
                        Dotplot_search_filter = paste0("H:/PCD_database/result3/",
                                                       cancer_names,
                                                       "/result_cpdb/Dotplot_search_filter.pdf"),
                        chord_plot_25 = paste0("H:/PCD_database/result3/",
                                              cancer_names,
                                              "/result_cpdb/chord_plot_25.pdf"),
                        chord_plot_05 = paste0("H:/PCD_database/result3/",
                                               cancer_names,
                                               "/result_cpdb/chord_plot_05.pdf"),
                        
                        Dotplot_search_filter_png = paste0("H:/PCD_database/result3_database/",
                                                       cancer_names,
                                                       "/Dotplot_search_filter.png"),
                        chord_plot_25_png = paste0("H:/PCD_database/result3_database/",
                                               cancer_names,
                                               "/chord_plot_25.png"),
                        chord_plot_05_png = paste0("H:/PCD_database/result3_database/",
                                               cancer_names,"/chord_plot_05.png"))
for(i in 1:30){
  dir.create(paste0("H:/PCD_database/result3_database/",cancer_names[i]))
  p = image_read_pdf(path = plot_file$Dotplot_search_filter[i])
  image_write(p,path =plot_file$Dotplot_search_filter_png[i])
  p = image_read_pdf(path = plot_file$chord_plot_25[i])
  image_write(p,path =plot_file$chord_plot_25_png[i])
  p = image_read_pdf(path = plot_file$chord_plot_05[i])
  image_write(p,path =plot_file$chord_plot_05_png[i])
}


#######################每个癌症的CPDB的互作基因结果

search_result_25
data_files <- paste0("H:/PCD_database/result3/",cancer_names,"/result_cpdb/search_result_25.csv")
interaction_res <- list()
for(i in 1:30 ){
  temp <-  read.csv(data_files[i])
  temp$cancer <- cancer_names[i]
  interaction_res[[i]] <- temp
}
interaction_res <- do.call(rbind,interaction_res)
interaction_res <- interaction_res[,c(10,8,1,2,3,4,5,6,7,9)]
write.csv(interaction_res,file = "H:/PCD_database/result3_database/interaction.csv")


######################结果三的泛癌统计
rm(list= ls())
load("互作统计.RData")
######################画那个大热图
cancer_names <- list.files("H:/PCD_database/result3/")
######################整理
data_files <- paste0("H:/PCD_database/result3/",cancer_names,"/result_cpdb")
interaction_number_fill <- list() 
interaction_sum <- list()
for(i in 1:30 ){
  file_temp <- list.files(data_files[i])
  file_temp_1 <- paste0(data_files[i],"/",grep("pvalues",file_temp,value = T))
  pvals <- read.delim(file_temp_1,
                      header = T,check.names = F )
  colnames(pvals) <- gsub("nan","None",x=colnames(pvals))
  pvals <- pvals %>% dplyr::filter(directionality=="Ligand-Receptor")
  interaction_number <- list()
  for(j in 14:dim(pvals)[2]){
    temp1 <-colnames(pvals)[j]
    temp2 <- count(pvals[,j] < 0.05)
    temp3 <- str_split(temp1,"\\|")[[1]][1]
    temp4 <- str_split(temp1,"\\|")[[1]][2]
    temp5 <- data.frame(interaction = temp1,
                        interaction_A = temp3,
                        interaction_B = temp4,
                        count = temp2,
                        cancer = cancer_names[i])
    interaction_number[[temp1]] <- temp5
  }
  
  temp_6 <- do.call(rbind,interaction_number)
  rownames(temp_6) <- NULL
  temp_6 <- temp_6[,c(2,3,4)]
  colnames(temp_6) <- c("SOURCE", "TARGET", "COUNT")
  count_mat  <- reshape2::acast(SOURCE ~ TARGET, data = temp_6, value.var = "COUNT")
  dcm <- diag(count_mat)
  count_mat <- count_mat + t(count_mat)
  diag(count_mat) <- dcm
  all_sum <- rowSums(count_mat)
  all_sum <- data.frame(all_sum)
  all_sum$cacner  <- cancer_names[i]
  all_sum$death <- rownames(all_sum)
  temp_6$cancer <- cancer_names[i]
  interaction_number_fill[[i]] <- temp_6
  interaction_sum[[i]] <-  all_sum
  
}

interaction_number_fill_data <- do.call(rbind,interaction_number_fill)

interaction_number_fill_data$he <- paste0(interaction_number_fill_data$SOURCE,"-",
                                          interaction_number_fill_data$TARGET)
  
interaction_sum_fill  <- do.call(rbind,interaction_sum)

#####################对数据框进行处理
interaction_number_fill_matrix <- matrix(0,ncol = 226,nrow = 30)

interaction_sum_fill_matrix <- matrix(0,ncol = 16,nrow = 30)

rownames(interaction_number_fill_matrix) <- cancer_names
colnames(interaction_number_fill_matrix) <- unique(interaction_number_fill_data$he)

rownames(interaction_sum_fill_matrix) <- cancer_names
colnames(interaction_sum_fill_matrix) <- unique(interaction_sum_fill$death)

for(i in 1:30){
  print(cancer_names[i])
  temp <- interaction_number_fill_data %>% filter(cancer == cancer_names[i])
  for(j in 1:dim(temp)[1]){
    interaction_number_fill_matrix[cancer_names[i],
                                   temp[j,5] ] <- temp[j,3]
  }
}

for(i in 1:30){
  print(cancer_names[i])
  temp <- interaction_sum_fill %>% filter(cacner == cancer_names[i])
  for(j in 1:dim(temp)[1]){
    interaction_sum_fill_matrix[cancer_names[i],
                                temp[j,3]] <-  temp[j,1]
  }
}

save(interaction_sum_fill_matrix,
     interaction_number_fill_matrix,file = "H:/PCD_database/结果三统计/统计.RData" )

load("H:/PCD_database/结果三统计/统计.RData")
#################画热图
library(pheatmap)
interaction_number_fill_matrix_1 <- scale(interaction_number_fill_matrix)
#interaction_number_fill_matrix_1 <- interaction_number_fill_matrix_1[,sort(colnames(interaction_number_fill_matrix_1))]
interaction_number_fill_matrix_1[is.na(interaction_number_fill_matrix_1)] <- 0

interaction_number_fill_matrix <- interaction_number_fill_matrix[,sort(colnames(interaction_number_fill_matrix ))]

#############前 113
dev.new()

p1 <- pheatmap(t(interaction_number_fill_matrix)[1:113,],cluster_rows = T,cluster_cols = F,
         color = colorRampPalette(c("#439551","#439551","white","firebrick3","firebrick3"))(100)
         ,display_numbers = F,scale="column",
         show_rownames = T,
         show_colnames = T,
         border= T,border_color = "#f5ece7")
############后 113
p2 <- pheatmap(t(interaction_number_fill_matrix)[114:226,],cluster_rows = T,cluster_cols = F,
         color = colorRampPalette(c("#439551","#439551","white","firebrick3","firebrick3"))(100)
         ,display_numbers = F,scale="column",
         show_rownames = T,
         show_colnames = T,
         border= T,border_color = "#f5ece7")

########## 每个死亡的互作强度
pheatmap(t(interaction_sum_fill_matrix),cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("#439551","#439551","white","firebrick3","firebrick3"))(100)
         ,display_numbers = F,scale="column",
         show_rownames = T,
         show_colnames = T,
         border= T,border_color = "#f5ece7")
max(interaction_sum_fill_matrix)

##########重要互作强度
pheatmap(t(interaction_number_fill_matrix)[unique(paste0(interaction_data$A,"-",
                                                  interaction_data$B)),],
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("#439551","#439551","white","firebrick3","firebrick3"))(100)
         ,display_numbers = F,scale="column",
         show_rownames = T,
         show_colnames = T,
         border= T,border_color = "#f5ece7")

##########重要互作统计
interaction_data$cancer[c(9,10)] <- "BRCA_TNBC" 
interaction_data_matrix <- matrix(0,nrow = 30,
                                  ncol = length(unique(interaction_data$interaction)))
colnames(interaction_data_matrix) <- unique(interaction_data$interaction)
rownames(interaction_data_matrix) <- cancer_names

for(i in 1:30){
  print(cancer_names[i])
  temp <- interaction_data %>% filter(cancer  == cancer_names[i])
  temp <- as.data.frame(temp)
  for(j in 1:dim(temp)[1]){
    interaction_data_matrix[cancer_names[i],
                                temp[j,2]] <- 1
  }
}


pheatmap(t(interaction_number_fill_matrix_1)[unique(paste0(interaction_data$A,"-",
                                                           interaction_data$B)),],
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("navy","white","firebrick3"))(100)
         ,display_numbers = F,scale="column",
         show_rownames = T,
         show_colnames = T,
         border= F)

pheatmap(t(interaction_data_matrix),
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("navy","white","firebrick3"))(100)
         ,display_numbers = F,scale="column",
         show_rownames = T,
         show_colnames = T,
         border= T)

###################读取significant_means.txt并得到泛癌的数据
rm(list=ls())
cancer_names <- list.files("H:/PCD_database/result3/")
load("互作统计.RData")
interaction_data$cancer[c(9,10)] <- "BRCA_TNBC"
##################
interaction_data$interaction_names_a  <- paste0(interaction_data$A,"|",interaction_data$B)
data_list <- list()
interaction_names_a <- unique(interaction_data$interaction_names_a)
for(i in interaction_names_a){
  temp <- which(interaction_data$interaction_names_a == i)
  temp1 <- unique(interaction_data$cancer[temp])
  temp_list <- list()
  for(j in 1:length(temp1)){
    temp_file <- paste0("H:/PCD_database/result3/",temp1[j],"/result_cpdb")
    temp_file <- grep("statistical_analysis_significant_means",
                      list.files(temp_file),value = T)
    temp_file <- paste0("H:/PCD_database/result3/",temp1[j],"/result_cpdb/",temp_file)
    interaction_means <- read.table(temp_file,header = T,fill = T,sep = "\t")
    colnames(interaction_means) <- gsub(pattern = "\\.",replacement = "|",colnames(interaction_means))
    colnames(interaction_means) <- gsub(pattern = "nan",replacement = "None",colnames(interaction_means))
    temp_111 <- which(colnames(interaction_means)==i)
    interaction_means <- interaction_means[,c(1:14,temp_111)]
    colnames(interaction_means)[15] <- "Means"
    interaction_means$interaction <- i
    interaction_means$cancer <- temp1[j]
    temp_list[[j]] <- interaction_means
  }
  temp_22222 <- do.call(rbind,temp_list)
  data_list[[i]] <- temp_22222
}

data_interaction <- do.call(rbind,data_list)
data_interaction <- data_interaction[-which(is.na(data_interaction$Means)),]



write.csv(data_interaction,"H:/PCD_database/结果三统计/interaction_patterns_ligand_receptor.csv")
####
temp <- read.csv("H:/PCD_database/结果三统计/interaction_patterns_ligand_receptor.csv")
unique(temp$interaction)


temp1 <- temp[which(temp$interaction=="Pyroptosis-Immunogenic_cell_death"),]

temp1 <- temp1[,c(4,15,17,19)]
temp1 <- temp1[-which(temp1$classification==""),]

temp2 <- temp1[,c(2,3,4)]
ccc <- aggregate(temp2$Means, by = list(temp2$classification,temp2$cancer), sum)
ccc <- as.data.frame(ccc)


ccc <- ccc[which(ccc$Group.1%in%c(
  "Signaling by Galectin"  ,                                        
  "Signaling by Lipoxin/Leukotriene" ,                              
  "Signaling by HLA"          ,                                     
  "Signaling by Selectin"       ,                                   
  "Adhesion by ICAM"    ,
  "Signaling by Amyloid-beta precursor protein",
  "Signaling by Chemokines"                ,                        
  "Signaling by Cholesterol/Desmosterol"       ,                    
  "Signaling by Epidermal growth factor"     ,                      
  "Signaling by Growth arrest"              ,                       
  "Signaling by Insulin-like growth factor"     ,                   
  "Signaling by Notch"                        ,                     
  "Signaling by Tumor necrosis factor"    ,                         
  "Signaling by Vascular endothelial growth factor"     ,           
  "Signaling by WNT"                )),]
unique(ccc$Group.1)
library(ggplot2)
ggplot(ccc)+geom_bar(aes(x=Group.2,y=x,fill= Group.1),stat = "identity")+
  scale_fill_manual(values = c("#FFB6C1", "#D3D3D3", "#8B4513", "#87CEFA", "#FFD700", 
                               "#98FB98", "#FF6347", "#20B2AA", "#B0E0E6", "#C71585", 
                               "#DAA520", "#FF4500", "#D2691E", "#FF1493", "#32CD32", 
                               "#8A2BE2", "#F0E68C", "#FF8C00", "#C0C0C0", "#F08080", 
                               "#7FFF00", "#40E0D0", "#ADFF2F", "#9932CC", "#6495ED", 
                               "#FF7F50", "#1E90FF", "#B22222", "#FFD700", "#F4A300"))+
  theme_classic()

    
ggplot(data = ccc) + geom_bar(mapping = aes(x=Group.2,y=x,fill= Group.1), stat = 'identity') +
  coord_polar(theta = "x")+
  scale_fill_manual(values = c("#FFB6C1", "#D3D3D3", "#8B4513", "#87CEFA", "#FFD700", 
                               "#98FB98", "#FF6347", "#20B2AA", "#B0E0E6", "#C71585", 
                               "#DAA520", "#FF4500", "#D2691E", "#FF1493", "#32CD32", 
                               "#8A2BE2", "#F0E68C", "#FF8C00", "#C0C0C0", "#F08080", 
                               "#7FFF00", "#40E0D0", "#ADFF2F", "#9932CC", "#6495ED", 
                               "#FF7F50", "#1E90FF", "#B22222", "#FFD700", "#F4A300"))+
  theme_classic()






"Signaling by Galectin"                                          
"Signaling by Lipoxin/Leukotriene"                               
"Signaling by HLA"                                               
"Signaling by Selectin"                                          
"Adhesion by ICAM"    
"Signaling by Amyloid-beta precursor protein"
"Signaling by Chemokines"                                        
"Signaling by Cholesterol/Desmosterol"                           
"Signaling by Epidermal growth factor"                           
"Signaling by Growth arrest"                                     
"Signaling by Insulin-like growth factor"                        
"Signaling by Notch"                                             
"Signaling by Tumor necrosis factor"                             
"Signaling by Vascular endothelial growth factor"                
"Signaling by WNT"                          
