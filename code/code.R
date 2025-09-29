library(Seurat)

meta_add <- read.csv("CART_All.integrated_metadata.csv",row.names=1)

seurat_obj <- readRDS("GSE262072_CART_Processed_object.rds")

seurat_obj@meta.data$Patient_ID <- substr(seurat_obj@meta.data$Sample_ID,1,4)

seurat_obj@meta.data$integrated_snn_res.0.4 <- meta_add[,"integrated_snn_res.0.4"]
seurat_obj@meta.data$seurat_clusters <- meta_add[,"seurat_clusters"]
seurat_obj@meta.data$Cell_type <- meta_add[,"Cell_type"]
seurat_obj@meta.data$condition <- substr(meta_add[,"patient"],1,2)



cell_cluster_annotation <- meta_add[,c("seurat_clusters","Cell_type")]
cell_cluster_annotation <- cell_cluster_annotation[!duplicated(cell_cluster_annotation[,1]),]
cell_cluster_annotation <- cell_cluster_annotation[order(cell_cluster_annotation[,1]),]
rownames(cell_cluster_annotation) <- NULL



patient_cell_count <- data.frame(table(seurat_obj@meta.data$Patient_ID))
colnames(patient_cell_count) <- c("patient","count")
write.csv(patient_cell_count,"patient_cell_count.csv",quote=F)



counts <- seurat_obj@assays$RNA@counts

patients <- as.character(patient_cell_count[,1])
for(patient in patients)
{
	print(patient)
	pseudobulk <- data.frame()

	counts_sub <- counts[,seurat_obj@meta.data$Patient_ID==patient]
	meta_sub <- seurat_obj@meta.data[seurat_obj@meta.data$Patient_ID==patient,]
	
	# extract count matrix
	expr <- as.matrix(counts_sub)
	colnames(expr) <- colnames(counts_sub)
	rownames(expr) <- rownames(counts_sub)
	
	# normalize to TPM
	expr <- t(t(expr)*1e5/colSums(expr))
	
	# transform to log space
	expr <- log2(expr + 1)
	
	celltypes <- unique(meta_sub[,"Cell_type"])
	for(celltype in celltypes)
	{
		expr_celltype <- expr[,meta_sub[,"Cell_type"]==celltype,drop=F]
		pseudobulk[rownames(expr_celltype),paste0(celltype)] <- apply(expr_celltype,1,mean)
	}
	
	conditions <- unique(meta_sub[,"condition"])
	for(condition in conditions)
	{
		expr_condition <- expr[,meta_sub[,"condition"]==condition,drop=F]
		pseudobulk[rownames(expr_condition),paste0(condition)] <- apply(expr_condition,1,mean)
	}
	
	write.csv(pseudobulk,paste0(patient,".csv"),quote=F)
}

