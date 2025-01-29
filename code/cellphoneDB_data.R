# 循环result1的癌症数据输出为cellphoneDB的准备文件###########
library(Seurat)
library(SeuratDisk)
library(readr)
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
devtools::install_local("~/R/x86_64-pc-linux-gnu-library/4.2/ktplots-master.zip")
library(ktplots)
library(SingleCellExperiment)

# 设定result1为base directory
setwd("~/ouqi/result1/")
dir.create("../result3/")
for (i in 1:length(list.files())) {
  # 创建counts_normalized h5ad
  load(paste0(list.files()[i],"/","PCD_ANN.RData"))
  #anndata <- new_data@assays$RNA@data %>% as.data.frame()
  dir.create(paste0("../result3/",list.files()[i]))
  # R
  # take raw data and normalise it
  count_raw <- new_data@assays$RNA@counts
  count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
  write.table(count_norm, paste0("../result3/",list.files()[i],"/CellphoneDB_count.txt"), sep='\t', quote=F)
  # 创建细胞metadata
  metadata <- new_data@active.ident %>% as.data.frame()
  colnames(metadata) <- "Death_state"
  metadata$barcode_sample <- rownames(metadata)
  metadata <- metadata %>% dplyr::select(2,1)
  write_tsv(metadata,paste0("../result3/",list.files()[i],"/metadata.tsv"))
  rm(list = ls())
  gc()
}
dir.create("CPDB")
# 开始cellphoneDB#####
# 进入linux系统

Linux+python代码在另一处




##########细胞互作##########
#######Back to R ###########
rm(list = ls())
setwd("~/ouqi/result3/")
cancer_name <- "PRAD"
# 读入cellphoneDB跑出来的文件
tmp <- list.files(paste0(cancer_name,"/result_cpdb/"))
pvals <- read.delim(paste0(cancer_name,"/result_cpdb/",tmp[grep(pattern = "pvalues",tmp)]),header = T,check.names = F)
colnames(pvals) <- gsub("nan","None",x=colnames(pvals))
means <- read.delim(paste0(cancer_name,"/result_cpdb/",tmp[grep(pattern = "analysis_means",tmp)]),header = T,check.names = F)
colnames(means) <- gsub("nan","None",x=colnames(means))
#interaction_scores <- read.delim(paste0(cancer_name,"/result_cpdb/",tmp[grep(pattern = "interaction_scores",tmp)]),header = T,check.names = F)

# 只要Ligand-Receptor的
pvals <- pvals %>% dplyr::filter(directionality=="Ligand-Receptor")
means <- means %>% dplyr::filter(directionality=="Ligand-Receptor")
#interaction_scores <- interaction_scores %>% dplyr::filter(directionality=="Ligand-Receptor")

# 导入seurat对象
load(paste0("../result1/",cancer_name,"/PCD_ANN.RData"))
new_data$Death_state <- new_data@active.ident %>% as.character()


# 先画热图
heatmap <- plot_cpdb_heatmap(pvals = pvals, cellheight = 30, cellwidth = 30,main = "Cell-Cell Interactions between Death States")
pdf(paste0(cancer_name,"/result_cpdb/","heatmap.pdf"),width = 11.69,height = 8.27)
print(heatmap)
dev.off()



#参数：根据interaction score最小得分大于50 （query_minimum_score = 50）
# minimum score that an interaction must have to be filtered.
search_result <- read.table(paste0(cancer_name,"/result_cpdb/",tmp[grep(pattern = "search_results",tmp)]), sep="\t",header = T)
search_result <- search_result[,-1]
# search全部结果
means_search <- means[which(means$interacting_pair%in%search_result$interacting_pair),]
pvals_search <- pvals[which(pvals$interacting_pair%in%search_result$interacting_pair),]
# search结果前25%
quantile_score <- quantile(search_result$significant_mean)
search_result_25 <- search_result[which(search_result$significant_mean>quantile_score['75%']),]
means_25 <- means[which(means$interacting_pair%in%search_result_25$interacting_pair),]
pvals_25 <- pvals[which(pvals$interacting_pair%in%search_result_25$interacting_pair),]



# 先画热图
#heatmap <- plot_cpdb_heatmap(pvals = pvals, cellheight = 30, cellwidth = 30,main = "Cell-Cell Interactions between Death States")
#pdf(paste0(cancer_name,"/result_cpdb/","heatmap.pdf"),width = 11.69,height = 8.27)
#print(heatmap)
#dev.off()


# 总图--没做任何处理的：cpdb_statistical_analysis_method.call
Dotplot1 <- plot_cpdb(
  scdata=new_data,
  cell_type1=".",
  cell_type2=".",  # this means all cell-types
  celltype_key="Death_state",
  #idents = 'Death_state',
  means=means,
  pvals=pvals,
  #genes=c("PTPRC", "TNFSF13"),
  #title="Ligand-Receptor Pairs between Death States",
  max_size = 2,
  #interaction_scores =  ,
  keep_significant_only = F
)
pdf(paste0(cancer_name,"/result_cpdb/","Dotplot_raw.pdf"),width = 30,height = 200)
print(Dotplot1)
dev.off()

# 总图--search result：search_utils.search_analysis_results
# 参数：根据interaction score最小得分大于50 （query_minimum_score = 50）
Dotplot2 <- plot_cpdb(
  scdata=new_data,
  cell_type1=".",
  cell_type2=".",  # this means all cell-types
  celltype_key="Death_state",
  #idents = 'Death_state',
  means=means_search,
  pvals=pvals_search,
  #genes=c("PTPRC", "TNFSF13"),
  #title="Ligand-Receptor Pairs between Death States",
  max_size = 2,
  #interaction_scores =  ,
  keep_significant_only = F
)
pdf(paste0(cancer_name,"/result_cpdb/","Dotplot_search_all.pdf"),width = 20,height = 25)
print(Dotplot2)
dev.off()


# 根据interaction score最小得分大于50 （query_minimum_score = 50）
# 筛选互作强度高的配受体对：significant means排序前25%的配受体对保留
Dotplot3 <- plot_cpdb(
  scdata=new_data,
  cell_type1=".",
  cell_type2=".",  # this means all cell-types
  celltype_key="Death_state",
  #idents = 'Death_state',
  highlight_size=1,  #显著性红圈大小
  means=means_25,
  pvals=pvals_25,
  #genes=c("PTPRC", "TNFSF13"),
  #title="Ligand-Receptor Pairs between Death States",
  max_size = 5,  #点的整体大小
  #interaction_scores =  ,
  keep_significant_only = F
)
pdf(paste0(cancer_name,"/result_cpdb/","Dotplot_search_filter.pdf"),width = 15,height = 10)
print(Dotplot3)
dev.off()

# 开始制作弦图
# 根据interaction score最小得分大于50 （query_minimum_score = 50）
# 筛选互作强度高的配受体对：significant means排序前25%的配受体对保留
decon <- read.delim(paste0(cancer_name,"/result_cpdb/",tmp[grep(pattern = "analysis_deconvoluted_05",tmp)]), check.names = FALSE)
decon2 <- decon %>% dplyr::filter(id_cp_interaction%in%unique(pvals_25$id_cp_interaction))
# 转为SCE对象
new_data_SCE <- as.SingleCellExperiment(new_data)

# 弦图--前25%配受体对，主要展示细胞
pdf(paste0(cancer_name,"/result_cpdb/","chord_plot_25.pdf"),width = 20,height = 10)
plot_cpdb3(
  scdata = new_data_SCE,
  cell_type1 = ".",
  cell_type2 = ".",
  celltype_key = "Death_state", # column name where the cell ids are located in the metadata
  means = means_25,
  pvals = pvals_25,
  deconvoluted = decon2 # new options from here on specific to plot_cpdb3
)
dev.off()



# 弦图
#chord plot
# 重新计算前interaction score 前5%的配受体对
quantile_score2 <- quantile(search_result$significant_mean,0.95)
search_result_05 <- search_result[which(search_result$significant_mean>quantile_score2),]
means_05 <- means[which(means$interacting_pair%in%search_result_05$interacting_pair),]
pvals_05 <- pvals[which(pvals$interacting_pair%in%search_result_05$interacting_pair),]
decon3 <- decon %>% dplyr::filter(id_cp_interaction%in%unique(pvals_05$id_cp_interaction))

interaction_annotation_used <- read.table("../CPDB/interaction_annotation_used.txt",header = T)
tmp2 <- means_05 %>% dplyr::mutate(anno=case_when(
  !(interacting_pair%in%interaction_annotation_used$interaction) ~ "Unknown",
  TRUE ~ "Unknown"
)) %>% dplyr::select(interacting_pair,anno)
# 制作配受体基因的细胞交互弦图
chord_plot <- plot_cpdb2(
  scdata = new_data_SCE,
  cell_type1 = ".",
  cell_type2 = ".",
  celltype_key = "Death_state", # column name where the cell ids are located in the metadata
  means = means_05,
  pvals = pvals_05,
  deconvoluted = decon3,# new options from here on specific to plot_cpdb2
  interaction_grouping = tmp2,
  edge_group_colors =c("Activating" = "#e15759",
                       "Inhibitory" = "#4e79a7",
                       "Chemotaxis" = "#59a14f",
                       "Intracellular trafficking" = "#9c755f",
                       "Unknown" = "#c08eaf"),
  node_group_colors = c("Alkaliptosis"="#41AB5D",
                        "Anoikis" = "#A35D34",
                        "Autophagy"="#ABCADE","Cuproptosis"="#8DD3C7",
                        "Entotie_cell_death"="#B3DE69",
                        "Extrinsic_apoptosis"="#3B77AC",
                        "Ferroptosis"="#BEBADA",
                        "Immunogenic_cell_death"="#C6DBEF",
                        "Intrinsic_apoptosis"="#80B1D3",
                        'Lysosome_dependent_cell_death'="#c06f98",
                        "Necroptosis"="#FB8072",
                        "Netotic_cell_death"="#EF9BB5",
                        "Oxeiptosis"="#F3B161",
                        "Parthanatos"="#B7DA90","Pyroptosis" = "#EB9B98",
                        "Extrinsic_Intrinsic_apoptosis"="#CF352B",
                        "None"="#8A959B"),
  plot_score_as_thickness = F,
  standard_scale = F
)
pdf(paste0(cancer_name,"/result_cpdb/","chord_plot_05.pdf"),width = 15,height = 10)
print(chord_plot)
dev.off()

# 最后存表search_result_25.csv!!!!!!!!!!!!!!
write.table(search_result_25,paste0(cancer_name,"/result_cpdb/","search_result_25.csv"),sep = ",",quote = F,row.names = F)





