scpUpload_newMeta <- function(gcdata, file.prefix, metadata=NULL, metatype = NULL, clustercols, clustertype, clusterfileOnly=FALSE) {
  
  #' This function generates the files needed for upload to the SCP. It was updated 03/03/2021 to comply with new SCP metadata requirements. 
  #' This function takes a long time to run because the generation of an .mtx file is lengthy. Please be patient. 
  
  #' @param gcdata: Seurat object 
  #' @param file.prefix: string indicating prefix you'd like to attach to your SCP files
  #' @param metadata: data frame containing both the metadata required by the SCP and the additional metadata you would like to include. 
  #' The following metadata columns MUST exist for upload to succeed: 
  #'    biosample_id: sample ID 
  #'    donor_id: donor ID 
  #'    species: species ID (e.g., NCBITaxon_10090; ignore any mutations for now)
  #'    species__ontology_label: label (e.g., Mus musculus)
  #'    sex: sex of each donor 
  #'    average_intensity: example of an additional column of our choosing?? 
  #'    disease: diease label from MONDO or PATO (e.g., MONDO_0005062, no disease is PATO_0000461)
  #'    diease__ontology_label: label for disease from MONDO or PATO (e.g., lymphoma, no disease is normal)
  #'    organ: ID for organ from Uberon (e.g., UBERON_0000029)
  #'    organ__ontology_label: label for organ (e.g., lymph node) 
  #'    library_preparation_protocol: ID for library prep (e.g., EFO:0009922)
  #'    library_preparation_protocol__ontology_label: label for protocol (e.g., 10x 3' v3 sequencing)
  #' The following are optional metadata that improve searchability of SCP portal studies
  #'    cell_type: OPTIONAL, ontology ID for cell type from EBI (e.g., CL_0000084)
  #'    cell_type__ontology_label: OPTIONAL cell type label (e.g., T cell)
  #' @param metatype: vector of strings indicating whether the metadata column is numeric or a "group" variable; the first entry MUST say "TYPE"
  #' For example: metatype <- c("TYPE", rep("group", 11), rep("numeric", 4), rep("group", 1))
  #' @param clustercols: a vector of strings indicating column names from your Seurat object metadata that you would like to have in your clustering file; 
  #' You can have multiple cluster files, one per subset of cells from main gene expression matrix
  #' For example, you might one one cluster file for all cells, and then one for subclustered T cells 
  #' @param clustertype: a vector of strings indicating whether the cluster column is numeric or a "group" variable 
  #' @param clusterfileOnly: FALSE if you'd like to also generate metadata and gene expression matrix files;
  #'                         TRUE if you only need to generate additional clustering files (e.g., for sub-clustering analyses)
  #'                         If clusterfileOnly=TRUE, gcdata (your Seurat object) should only contain the subset of cells you want associated 
  #'                         with this additionalcluster file, and the UMAP coordinates in that Seurat object should match the sub-analysis 
  #' @return Five files are returned if clusterfileOnly=FALSE: an _expression_data.mtx file containing the expression data, an _geneNames.csv file with 
  #' gene names, a _cellNames.csv with cell names, a _metaData.txt file containing metadata, and a _clusterfile.txt containing the cluster info 
  #' associated with the UMAP coordinates. If clusterfileOnly=TRUE, only a single file, the _clusterfile.txt file, will be returned.                          
  
  ######## write expression matrix as .mtx file ###############
  # After running this step, you need to gene sort using the SortSparseMatrix.py file provided by SCP 
  if (!clusterfileOnly){
    # library(Matrix)
    # writeMM(gcdata[['RNA']]@data, file=paste0(file.prefix, '_expression_data.mtx'))
    # 
    # ###### create gene and barcode name files, to be "bundled" with the .mtx file ###############
    # write.table(rownames(gcdata[['RNA']]@data), file=paste0(file.prefix, '_geneNames.csv'), row.names = FALSE, col.names=FALSE, sep=",", quote=FALSE)
    # write.table(colnames(gcdata[['RNA']]@data),file=paste0(file.prefix, '_cellNames.csv'), row.names = FALSE, col.names = FALSE, sep=",", quote=FALSE)

    ######## write meta data ############
    metadata$NAME <- rownames(metadata) # format for SCP
    metadata <- metadata[,c('NAME', setdiff(colnames(metadata), 'NAME'))] # format for SCP
    metadata[] <- lapply(metadata, as.matrix) # make sure none of the data frame columns are factors
    metadata[1,] <- metatype # set the type for each metadata column (numeric or group) for SCP visualization
    write.table(metadata, file=paste0(file.prefix, '_metaData.txt'), sep="\t", row.names = FALSE, quote=FALSE)
  }
  
  ####### create cluster file ############
  clusterfile <- as.data.frame(gcdata@reductions$umap@cell.embeddings) #extract umap coordinates 
  cluster.metadata <- gcdata@meta.data[, clustercols]
  if (length(clustercols)==1){cluster.metadata <- as.matrix(cluster.metadata)}
  if (length(clustercols)>1){lapply(cluster.metadata, as.matrix)}
  colnames(cluster.metadata) <- clustercols 
  clusterfile <- cbind(clusterfile, cluster.metadata, stringsAsFactors=FALSE)
  clusterfile$NAME <- rownames(clusterfile) # format for SCP 
  clusterfile <- clusterfile[,c('NAME', setdiff(colnames(clusterfile), 'NAME'))] # format for SCP 
  clusterfile[1,] <- cbind("TYPE", "numeric", "numeric", t(as.matrix(clustertype))) # format for SCP / set column types for the coordinates 
  colnames(clusterfile)[1:3] <- c("NAME", "X", "Y") # format for SCP 
  write.table(clusterfile, file=paste0(file.prefix, '_clusterfile.txt'), sep="\t", row.names = FALSE, quote=FALSE)
  
}

scpUpload <- function(gcdata, file.prefix, metacols=NULL, metatype=NULL, clustercols, clustertype, clusterfileOnly=FALSE) {
  
  ######## write expression matrix as .mtx file ###############
  # After running this step, you need to gene sort using the SortSparseMatrix.py file provided by SCP 
  if (!clusterfileOnly){
    writeMM(gcdata[['RNA']]@data, file=paste0(file.prefix, '_expression_data.mtx'))

    ###### create gene and barcode name files, to be "bundled" with the .mtx file ###############
    write.table(rownames(gcdata[['RNA']]@data), file=paste0(file.prefix, '_geneNames.csv'), row.names = FALSE, col.names=FALSE, sep=",", quote=FALSE)
    write.table(colnames(gcdata[['RNA']]@data),file=paste0(file.prefix, '_cellNames.csv'), row.names = FALSE, col.names = FALSE, sep=",", quote=FALSE)

    ######## write meta data ############
    metadata <- gcdata@meta.data[, metacols] # extract columns of interest
    metadata$NAME <- rownames(metadata) # format for SCP
    metadata <- metadata[,c('NAME', setdiff(colnames(metadata), 'NAME'))] # format for SCP
    metadata[] <- lapply(metadata, as.matrix) # make sure none of the data frame columns are factors
    metadata[1,] <- metatype # set the type for each metadata column (numeric or group) for SCP visualization
    write.table(metadata, file=paste0(file.prefix, '_metaData.txt'), sep="\t", row.names = FALSE, quote=FALSE)
  }
  
  ####### create cluster file ############
  clusterfile <- as.data.frame(gcdata@reductions$umap@cell.embeddings) #extract umap coordinates 
  cluster.metadata <- gcdata@meta.data[, clustercols]
  if (length(clustercols)==1){cluster.metadata <- as.matrix(cluster.metadata)}
  if (length(clustercols)>1){lapply(cluster.metadata, as.matrix)}
  colnames(cluster.metadata) <- clustercols 
  clusterfile <- cbind(clusterfile, cluster.metadata, stringsAsFactors=FALSE)
  clusterfile$NAME <- rownames(clusterfile) # format for SCP 
  clusterfile <- clusterfile[,c('NAME', setdiff(colnames(clusterfile), 'NAME'))] # format for SCP 
  clusterfile[1,] <- cbind("TYPE", "numeric", "numeric", as.matrix(clustertype)) # format for SCP / set column types for the coordinates 
  colnames(clusterfile)[1:3] <- c("NAME", "X", "Y") # format for SCP 
  write.table(clusterfile, file=paste0(file.prefix, '_clusterfile.txt'), sep="\t", row.names = FALSE, quote=FALSE)
  
}

scpUpload_singleUMAP <- function(gcdata, file.prefix, metacols, metatype) {
  
  ######## write expression matrix as .mtx file ###############
  # After running this step, you need to gene sort using the SortSparseMatrix.py file provided by SCP 
  # writeMM(gcdata[['RNA']]@data, file=paste0(file.prefix, '_expression_data.mtx'))
  
  ####### create gene and barcode name files, to be "bundled" with the .mtx file ###############
  write.table(rownames(gcdata[['RNA']]@data), file=paste0(file.prefix, '_geneNames.csv'), row.names = FALSE, col.names=FALSE, sep=",", quote=FALSE)
  write.table(colnames(gcdata[['RNA']]@data),file=paste0(file.prefix, '_cellNames.csv'), row.names = FALSE, col.names = FALSE, sep=",", quote=FALSE)
  
  ######## write meta data ############
  metadata <- gcdata@meta.data[, metacols] # extract columns of interest
  metadata$NAME <- rownames(metadata) # format for SCP 
  metadata <- metadata[,c('NAME', setdiff(colnames(metadata), 'NAME'))] # format for SCP 
  metadata[] <- lapply(metadata, as.matrix) # make sure none of the data frame columns are factors 
  metadata[1,] <- metatype # set the type for each metadata column (numeric or group) for SCP visualization 
  write.table(metadata, file=paste0(file.prefix, '_metaData.txt'), sep="\t", row.names = FALSE, quote=FALSE)
  
  ####### create cluster file ############
  clusterfile <- as.data.frame(gcdata@reductions$umap@cell.embeddings) # extract umap coordinates 
  clusterfile$NAME <- rownames(clusterfile) # format for SCP 
  clusterfile <- clusterfile[,c('NAME', setdiff(colnames(clusterfile), 'NAME'))] # format for SCP 
  clusterfile[1,] <- c("TYPE", rep("numeric", 2)) # format for SCP / set column types for the coordinates 
  colnames(clusterfile) <- c("NAME", "X", "Y") # format for SCP 
  write.table(clusterfile, file=paste0(file.prefix, '_clusterfile.txt'), sep="\t", row.names = FALSE, quote=FALSE)
  
}

PlotDoubletEmptyDrop <- function(df, file.plot, title, h=7, w=9, pt.size=0.1) {
  df.true <- df[grep("True|TRUE|1", df$metadata),]
  df.false <- df[grep("False|FALSE|0", df$metadata),]
  p1 <- ggplot() +
    geom_point(data = df.false, mapping = aes(UMAP_1, UMAP_2), color="gray50", size = .5) +
    geom_point(data = df.true, mapping = aes(UMAP_1, UMAP_2), color = "red", size = .5) +
    labs(color = title, title=title) +
    scale_colour_manual(values = c("gray50", "red")) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(axis.text = element_text(size=8, family="Helvetica"),
          legend.position = "right",
          legend.text = element_text(size=8, family="Helvetica"),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,0,-10,-10),
          legend.title = element_text(size=8, family="Helvetica"),
          legend.key.size = unit(.1, "cm"),
          axis.line.x = element_line(color = 'black', size = 0.75),
          axis.line.y = element_line(color = 'black', size = 0.75),
          axis.ticks.x = element_line(color = 'black', size = 0.75),
          axis.ticks.y = element_line(color = 'black', size = 0.75),
          axis.title.y = element_text(size=8, family="Helvetica", margin = margin(0, 1, 0, 0)),
          axis.title.x = element_text(size=8, family="Helvetica", margin = margin(1, 0, 0, 0)),
          aspect.ratio = 1, plot.margin = unit(c(.1, .3, .1, .1), "cm"),
          plot.background = element_blank(),
          plot.title = element_text(size=8, family="Helvetica")) +
    guides(color = guide_legend(override.aes = list(size=2)))
  save_plot(file.plot, p1, base_height = h, base_width = w)
}

PlotMetaData <- function(df, file.plot, title, h=7, w=9, pt.size=0.1) {
  p1 <- ggplot(data = df) +
    geom_point(mapping = aes(UMAP_1, UMAP_2, color = metadata), size = pt.size, shape=19) +
    labs(color = title, title="") +
    scale_colour_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(axis.text = element_text(size=8, family="Helvetica"),
          legend.position = "right",
          legend.text = element_text(size=8, family="Helvetica"),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,0,-10,-10),
          legend.title = element_text(size=8, family="Helvetica"),
          legend.key.size = unit(.1, "cm"),
          axis.line.x = element_line(color = 'black', size = 0.25),
          axis.line.y = element_line(color = 'black', size = 0.25),
          axis.ticks.x = element_line(color = 'black', size = 0.25),
          axis.ticks.y = element_line(color = 'black', size = 0.25),
          axis.title.y = element_text(size=8, family="Helvetica", margin = margin(0, 1, 0, 0)),
          axis.title.x = element_text(size=8, family="Helvetica", margin = margin(1, 0, 0, 0)),
          aspect.ratio = 1, plot.margin = unit(c(.1, .3, .1, .1), "cm"),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          #             panel.grid.minor = element_blank(),
          plot.title = element_text(size=8, family="Helvetica")) +
    guides(color = guide_legend(override.aes = list(size=2)))
  save_plot(file.plot, p1, base_height = h, base_width = w)
}


PlotMetaDataBar <- function(df, file.plot, title, h=2, w=5) {
  p <- ggplot(df) + 
    geom_bar(mapping = aes(x=clusters, fill = metadata), position="fill") +           
    scale_fill_manual(values = colors) + 
    theme(axis.text = element_text(size=6, family="Helvetica"),
          legend.position = "right",
          legend.text = element_text(size=6, family="Helvetica"),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,0,-10,-10),
          legend.title = element_text(size=6, family="Helvetica"),
          legend.key.size = unit(.2, "cm"),
          axis.line.x = element_line(color = 'black', size = 0.25),
          axis.line.y = element_line(color = 'black', size = 0.25),
          axis.ticks.x = element_line(color = 'black', size = 0.25),
          axis.ticks.y = element_line(color = 'black', size = 0.25),
          axis.title.y = element_text(size=6, family="Helvetica", margin = margin(0, 1, 0, 0)),
          axis.title.x = element_text(size=6, family="Helvetica", margin = margin(1, 0, 0, 0)),
          aspect.ratio = 0.5, plot.margin = unit(c(.1, .3, .1, .1), "cm"),
          plot.background = element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          plot.title = element_text(size=6, family="Helvetica")) +
    labs(fill = title, title="") 
  save_plot(file.plot, p, base_height = h, base_width = w)
  print(p)
  return(p)
}

library(data.table)
# source the dirichlet_regression function 
source(paste0(user.path, '/code/dirichlet_regression.R'))
metadata_pval <- function(metadata, group, idents, levels, string, p.file) {
  # make counts matrix
  counts <- do.call(cbind,tapply(paste0(metadata, '_', group), idents, table))
  # make covariate data.frame 
  covariates = data.frame(condition=gsub(string, '', rownames(counts)))
  covariates$condition <- factor(covariates$condition, levels = levels)
  # run the regression to determine frequencies 
  res = dirichlet_regression(counts, covariates, counts ~ condition)
  res$padj = p.adjust(res$pval, method="fdr")
  write.csv(rbind(res$pvals, res$padj), p.file)
}
