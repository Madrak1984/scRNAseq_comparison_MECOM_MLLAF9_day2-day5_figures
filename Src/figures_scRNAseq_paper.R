# command lines for Athimed
# must be added as code for the paper submission
# done By Jonathan Seguin, group of Prof. Schwaller, DBM and UKBB
# Wed Oct 22 13:19:25 2025


# preparation of the working environment ----------------------------------

## load the libraries ----------------------------------
library(edgeR)
library(tidyverse)
library(pheatmap)
library(writexl)
library(biomaRt)
library(SingleCellExperiment)



## load the functions ----------------------------------
# source("Src/DE_functions.R")


##### function to save figures for paper in png, svg and pdf formats
# filename, character value, name of the file to save without extension
# ggplot, ggplot object, contains the figure to save in the file
# dirPlot, character value (current working directory by default), indicate the path were the file must be saved
# A4, boolean value (FALSE by default), indicate if the figure must be saved in A4 landscape format
#####
saveFigures <- function(fileName, ggplot, dirPlot = getwd(), A4 = F){
  
  if(!dir.exists(dirPlot)){
    stop("dirPlot does not exist!")
  }
  
  # save the figures in png, svg and pdf format
  # save in the svg format
  if(A4){
    svg(filename = file.path(dirPlot, paste0(fileName, ".svg")), width = 11.7, height = 8.3)
  } else {
    svg(filename = file.path(dirPlot, paste0(fileName, ".svg")))
  }
  print(ggplot)
  dev.off()
  
  # save in the png format
  if(A4){
    png(filename = file.path(dirPlot, paste0(fileName, ".png")), width = 297, height = 210, units = "mm", res = 300)
  } else {
    png(filename = file.path(dirPlot, paste0(fileName, ".png")))
  }
  print(ggplot)
  dev.off()
  
  # save in the pdf format
  if(A4){
    pdf(file = file.path(dirPlot, paste0(fileName, ".pdf")), paper = "a4r", width = 11.7, height = 8.3)
  } else {
    pdf(file = file.path(dirPlot, paste0(fileName, ".pdf")))
  }
  print(ggplot)
  dev.off()
  
}


# function to create a directory
###### function to create a directory
# path, character value, indicate the path to create if the folders did not exist yet
##### return path
create_dir <- function(path){
  
  if(is.character(path)){
    
    if(is.vector(path)){
      
      if(!dir.exists(paths = path)) dir.create(path = path, recursive = T)
      
    } else stop("path must be a vector!")
    
  } else stop("path must be a character!")
  
  return(path)
}

# function to have updated cluster
getAnnotation <- function(cluster_name){
  
  # get the dataframe containing the annotation
  dt_annotation <- data.frame("Cluster" = paste0("clust", c("000", "001", "002", paste0("0", 1:9), 10:19)),
                              "Annotation" = c("ST-HSC_1", "LT-HSC", "ST-HSC_2", "MPP1_1", "MPP3_1",
                                               "MPP2_1", "MPP3_2", "MPP3_3", "MPP1_2", "MPP2_2", 
                                               "MPP3_4", "GMP_1", "GMP_2", "MPP2_3", "CMP_1", 
                                               "MEP", "T/NK cells pro", NA, "CMP_2", "CLP",
                                               "B cells pro", NA))
  
  # return the annotation 
  return(dt_annotation$Annotation[dt_annotation$Cluster == cluster_name])
  
  
}



######   function to draw the volcanoplot
#   dgeobject, DGELRT or dataframe object, contains the Differential expression analysis done by EdgeR
#   FC, numeric value, indicate a threshold for the FC (NULL by default, and not implemented yet)
#   FDR, numeric value, threshold to change the color according to the FDR (0.05 by default)
#   genes, vector of character, genes names (SYMBOL) to display on the volcanoplot (NULL by default, see addGenes parameters)
#   addGenes, boolean value, indicate if genes must be added (TRUE by default, add significant genes by default)
#   setAxis, boolean value, indicate if we want to set the FC axis to be symetrical (FALSE by defaults)
#   nbGenes, numeric value, indicate the number of significant genes to display (NULL by default)
#   yaxis, numeric values, vector object, must indicate the limit for the y-axis with 2 values (minimal and maximal), NULL by default
#   xaxis, numeric values, vector object, must indicate the limit for the x-axis with 2 values (minimal and maximal), NULL by default
#   displayPVal, boolean value, must indicate if we display the pval not the FDR for the y-axis (FALSE by default)
#   pval, numeric value, threshold to change the color according to the selected pvalue (NULL by default)
#####
getVolcanoplot <-  function(dgeobject, FC = NULL, FDR= 0.05, genes = NULL,  addGenes = T,
                            setAxis = F, nbGenes = NULL, yaxis = NULL, xaxis = NULL,
                            displayPVal = F, pval=NULL, colorFDR = T){
  
  #### check the parameters
  
  # check dgeobject
  if(class(dgeobject) == "DGELRT"){
    dge_DEG <- as.data.frame(topTags(dgeobject, n = nrow(dgeobject)))
  } else if(is.data.frame(dgeobject)){
    dge_DEG <- dgeobject
  } else stop("dgeobject is not a dataframe or a DGELRT object (see EdgeR package documentation for more details")
  
  # check FDR
  if(!is.numeric(FDR)) stop("FDR must be numeric")
  
  # check genes
  if(!is.logical(addGenes)) stop("addGenes must be boolean")
  if(addGenes & (!is.null(genes))){
    if(is.character(genes)){
      if(!is.vector(genes)) stop("genes must be a vector")
    } else stop("genes must contain character")
  }
  
  # check setAxis
  if(!is.logical(setAxis)) stop("setAxis must be boolean")
  
  # check the nbGenes
  if(addGenes & (!is.null(nbGenes))){
    if(!is.numeric(nbGenes)) stop("nbGenes must be numeric")
  }
  
  # check the y-axis
  if((!is.numeric(yaxis)) & (!is.null(yaxis))){
    stop("y-axis is not a nemurical values!")
  } else if((is.numeric(yaxis)) & (!is.null(yaxis))){
    if(length(yaxis) != 2){
      stop("yaxis must be a vector of length 2.")
    }
  }
  
  # add the colors
  dge_DEG$color <- "darkblue"
  if(displayPVal & (!is.null(pval)) & (!colorFDR)) {
    dge_DEG$color[dge_DEG$PValue <= pval] <- "red"
  } else {
    dge_DEG$color[dge_DEG$FDR <= FDR] <- "red"
  }
  
  
  # select genes to display
  if(addGenes){
    
    # select the genes
    if(is.null(genes)){
      
      if(displayPVal & (!is.null(pval))){
        all_genes <- dge_DEG$SYMBOL[dge_DEG$PValue <= pval]
      } else {
        all_genes <- dge_DEG$SYMBOL[dge_DEG$FDR <= FDR]
      }
    }
    
    # select genes according to the count
    if(!is.null(nbGenes)){
      
      if(displayPVal & (!is.null(pval))){
        
        # select the dataframe according to the FDR
        tmp_data <- dge_DEG[dge_DEG$PValue <= pval,]
        
        # select the ordered genes
        tmp_genes <- tmp_data$SYMBOL[order(tmp_data$PValue)] 
        
      } else {
        
        # select the dataframe according to the FDR
        tmp_data <- dge_DEG[dge_DEG$FDR <= FDR,]
        
        # select the ordered genes
        tmp_genes <- tmp_data$SYMBOL[order(tmp_data$FDR)]  
      }
      
      
      # remove genes already selected
      if(!is.null(genes)) tmp_genes <- tmp_genes[!(tmp_genes %in% genes)]
      
      # select the genes according to the nbGenes
      if(nbGenes < length(tmp_genes)){
        tmp_genes <- head(tmp_genes, n= nbGenes)
      }
      
      # concatenate the genes
      if(!is.null(genes)){
        genes <- c(genes, tmp_genes)
      } else genes <- tmp_genes
    } else if(is.null(genes)) genes <- all_genes
  }
  
  # create the ggplot
  if(displayPVal){
    volcano <- ggplot(dge_DEG, aes(x=logFC, y=-log10(PValue), color=color)) +
      theme_bw() + geom_point(alpha=0.5, color=dge_DEG$color)
  } else {
    volcano <- ggplot(dge_DEG, aes(x=logFC, y=-log10(FDR), color=color)) +
      theme_bw() + geom_point(alpha=0.5, color=dge_DEG$color)
  }
  
  
  # add genes
  if(addGenes){
    volcano <- volcano +
      geom_text_repel(data=subset(dge_DEG[dge_DEG$SYMBOL %in% genes,]), color="black", aes(label=SYMBOL), size =5)
  }
  
  
  # set the x-axis
  if(setAxis){
    
    # get the max values
    max_x <- ceiling(max(na.omit(abs(dge_DEG$logFC))))
    volcano <- volcano + xlim(-max_x, max_x)
  }
  #volcano <- volcano + ggtitle(title)
  
  # set the y-axis
  if(!is.null(yaxis)){
    volcano <- volcano + ylim(yaxis)
  }
  
  # set the x-axis
  if(!is.null(xaxis)){
    volcano <- volcano + xlim(xaxis)
  }
  
  
  # return the plot
  return(volcano)
}




##### get mouse gene annotation 
#
#####
getMouseGene <- function(){
  
  # load the ensembl gene information
  mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  
  # list all genes
  genes_mm <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"), 
                    mart = mart)
  
  return(genes_mm)
}


## set the R objects ----------------------------------

# load the edgeR object
edgeR_result <- readRDS("Datasets/edgeR_results.rds")

# create the output directory
outDir <- create_dir("Output/Plots")

# create the output directory for the DEA 
outDir_DEA <- create_dir("Output/DEA_Analysis/")



# figure 3 D & H ----------------------------------


## create a count table for dysregulated genes  ----------------------------------

# initiate the count table
dt_count_DEGs <- data.frame("comparisons" = names(edgeR_result),
                            "Up" = 0, 
                            "Down" = 0,
                            "total genes" = 0,
                            "cluster" =  sapply(names(edgeR_result), function(x)strsplit(x, split = "\\.")[[1]][1]),
                            "group" = sapply(names(edgeR_result), function(x)strsplit(x, split = "\\.")[[1]][2]),
                            "day" = sapply(names(edgeR_result), function(x)strsplit(x, split = "\\.")[[1]][3]),
                            check.names = F)



# count the DEG
for(comp in names(edgeR_result)){
  
  # print the starting message
  print(paste0(date(), "; start to count for the comparison ", comp))
  
  # select dysregulated genes
  dt_tmp <- as.data.frame(topTags(edgeR_result[[comp]], n = nrow(edgeR_result[[comp]])))
  
  # remove Riken cDNA and predicted genes
  dt_tmp <- dt_tmp[(!grepl(pattern = "predicted gene", dt_tmp$DESCRIPTION) & !grepl(pattern = "RIKEN", dt_tmp$DESCRIPTION)),]
  
  # update the count table
  dt_count_DEGs[dt_count_DEGs$comparisons %in% comp, "total genes"] = nrow(dt_tmp)
  dt_count_DEGs[dt_count_DEGs$comparisons %in% comp, "Up"] = filter(dt_tmp, FDR <= 0.05, logFC > 0) %>% summarise(count = n())
  dt_count_DEGs[dt_count_DEGs$comparisons %in% comp, "Down"] = filter(dt_tmp, FDR <= 0.05, logFC < 0) %>% summarise(count = n())
  
  # print the ending message
  print(paste0(date(), "; counting for the comparison ", comp, " done."))
  
}


## create the matrix table for the heatmap creation  ----------------------------------

# extract all the clusters name
clusters_name <- unique(dt_count_DEGs$cluster)

# create the column names for the matrix
comp_name <- paste(rep(paste0("day", c(2,5)), each = 2), c("Down", "Up"), sep = "_")

# create the matrix
mat_countDEG <- matrix(0, nrow = length(clusters_name),
                       ncol = length(comp_name),
                       dimnames = list(clusters_name, comp_name))

# update the matrix
mat_countDEG[filter(dt_count_DEGs, group == "fusion", day == "day2")%>% dplyr::select(cluster) %>% unlist(), "day2_Up"] =
  unlist(filter(dt_count_DEGs, group == "fusion", day == "day2") %>% arrange(cluster) %>% dplyr::select(Up))
mat_countDEG[filter(dt_count_DEGs, group == "fusion", day == "day5")%>% dplyr::select(cluster) %>% unlist(), "day5_Up"] =
  unlist(filter(dt_count_DEGs, group == "fusion", day == "day5") %>% arrange(cluster) %>% dplyr::select(Up))
mat_countDEG[filter(dt_count_DEGs, group == "fusion", day == "day2")%>% dplyr::select(cluster) %>% unlist(), "day2_Down"] =
  unlist(filter(dt_count_DEGs, group == "fusion", day == "day2") %>% arrange(cluster) %>% dplyr::select(Down))
mat_countDEG[filter(dt_count_DEGs, group == "fusion", day == "day5")%>% dplyr::select(cluster) %>% unlist(), "day5_Down"] =
  unlist(filter(dt_count_DEGs, group == "fusion", day == "day5") %>% arrange(cluster) %>% dplyr::select(Down))

# calculate the percentage
max_Down <- max(mat_countDEG[, grep(pattern = "Down", colnames(mat_countDEG))])
max_Up <- max(mat_countDEG[, grep(pattern = "Up", colnames(mat_countDEG))])
mat_countDEG[, grep(pattern = "Down", colnames(mat_countDEG))] = mat_countDEG[, grep(pattern = "Down", colnames(mat_countDEG))]*100/max_Down
mat_countDEG[, grep(pattern = "Up", colnames(mat_countDEG))] = mat_countDEG[, grep(pattern = "Up", colnames(mat_countDEG))]*100/max_Up


# reverse the values for the down DEG
mat_countDEG[, grep(pattern = "Down", colnames(mat_countDEG))] = 0 - mat_countDEG[, grep(pattern = "Down", colnames(mat_countDEG))]

# update the annotations
rownames(mat_countDEG) = sapply(rownames(mat_countDEG), getAnnotation)

# create the heatmap
heatmap <- pheatmap(mat = mat_countDEG, 
                    cluster_rows = F, 
                    cluster_cols = F, 
                    main = "# DEG for TPO vs PBS comparison (FDR <= 0.05)", 
                    angle_col = 0, 
                    color = colorRampPalette(colors = c("blue", "white", "darkred"))(100),
                    legend = T,
                    legend_breaks = c(-100, -50, 0, 50, 100),
                    legend_labels  = as.character(c(max_Down, round(max_Down/2), 0, round(max_Up/2), max_Up)))

# save the heatmap in files
saveFigures(fileName = "heatmap_countDEG_pval0.005", ggplot = heatmap, dirPlot = outDir)




# figure 3 E,F,I and J ----------------------------------


# create the volcanoplot for the cluster LT-HSC with specific genes at day5
# get the volcanoplot
volcan <- getVolcanoplot(dgeobject = edgeR_result$clust001.fusion.day5, displayPVal = T, 
                         addGenes = T, 
                         genes = c("Prmt6", "Pbx3", "Creb5", "Gm28438", "mt-Nd3", "Auts2", "Arhgap6"), yaxis = c(0,5))
# save the figure in the pdf file
saveFigures(fileName = "volcanoplot_LTHSC_fusion_day5_paper", ggplot = volcan, dirPlot = create_dir(file.path(outDir, "LT-HSC_fusion_day5")))


# create the volcanoplot for the cluster ST-HSC2 with specific genes at day2
# get the volcanoplot
volcan <- getVolcanoplot(dgeobject = edgeR_result$clust002.fusion.day2, displayPVal = T,
                         yaxis = c(0,10), pval = 0.005)
# save the figure in the pdf file
saveFigures(fileName = "volcanoplot_STHSC-2_fusion_day2_paper", ggplot = volcan, dirPlot = create_dir(file.path(outDir, "ST-HSC_2_fusion_day2")))


# create the volcanoplot for the cluster MPP1_2 at day2
# get the volcanoplot
volcan <- getVolcanoplot(dgeobject = edgeR_result$clust06.fusion.day2, displayPVal = T, 
                         pval = 0.005, addGenes = T, yaxis = c(0,7))
# save the figure in the pdf file
saveFigures(fileName = "volcanoplot_MPP1-2_fusion_day2_paper", ggplot = volcan, dirPlot = create_dir(file.path(outDir, "MPP1_2_fusion_day2")))


# create the volcanoplot for the cluster MPP2_3 at day5
# get the volcanoplot
volcan <- getVolcanoplot(dgeobject = edgeR_result$clust11.fusion.day5, displayPVal = T, 
                         pval = 0.005, addGenes = T, yaxis = c(0,11))
# save the figure in the pdf file
saveFigures(fileName = "volcanoplot_MPP2-3_fusion_day5_paper", ggplot = volcan, dirPlot = create_dir(file.path(outDir, "MPP2-3_fusion_day5")))





# figure 3 G (matplot_day2_fusion_FDR0.05 file)  and K (matplot_day5_fusion_FDR0.05 file) ----------------------------------


# load the camera result
camera_results <- readRDS(file = "Datasets/camera_results.rds")

# get the mouse genes
mouse_genes <- getMouseGene()

# init the list for camera results
list_cameras <- NULL

# save the different collection from Msig database for camera
db_cameras = c("H", "c2", "c3", "c4", "c5", "c6", "c7")

# start the loop
for(db in db_cameras){
  
  # create the filename
  filename = paste0("mouse_", db, "_v5p2.rdata")
  
  # create the file
  #file = file.path(getwd(), filename)
  
  # download the file
  if(!file.exists(filename)){
    print(paste0("download of ",  filename, " in progress..."))
    download.file(paste0("http://bioinf.wehi.edu.au/software/MSigDB/", filename), filename)
    print(paste0(filename, " downloaded."))
    load(filename)
    print(paste0(filename, " loaded."))
  }
  
}


# create the output directory
outputdir <- "Output/DEA_Analysis/Pathway_Analysis"
if(!dir.exists(outputdir)) dir.create(outputdir, recursive = T)


## select the pathways ----------------------------------

### genes selection ----------------------------------

# find row index of genes
id_Hoxa9 <- grep(pattern = "Hoxa9", mouse_genes$external_gene_name)
id_Myc <- grep(pattern = "Myc$", mouse_genes$external_gene_name)
id_Stat <- grep(pattern = "Stat", mouse_genes$external_gene_name)
id_Jak <- grep(pattern = "Jak[1-9]$", mouse_genes$external_gene_name)
id_Mapk <- grep(pattern = "Mapk", mouse_genes$external_gene_name)
ids <- c(id_Hoxa9, id_Myc, id_Jak, id_Stat, id_Mapk)

# selection of genes
selected_genes <- mouse_genes[ids,]


## pathway selection ----------------------------------

# create vector to save the pathway of interest
selected_pathways <- NULL

# selection of pathway which contains the genes
for(pathway in names(Mm.H)){
  
  if(any(Mm.H[[pathway]] %in% selected_genes$entrezgene_id)){
    print(paste0(pathway, " contains selected genes."))
    selected_pathways <- c(selected_pathways, pathway)
  }
}


# add oxidative pathway and TAKEDA
selected_pathways <- c(selected_pathways, "HALLMARK_OXIDATIVE_PHOSPHORYLATION", 
                       "TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_8D_UP",
                       "TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_16D_UP")


# create the heatmap ----------------------------------

# loop to analyze according to the day
for(day in c("day2", "day5")){
  
  # select the cluster of interest, ex: day 2
  results <- names(camera_results$h.CATEGORY.v7.5.1.ensembl.mmu)
  selected_cluster <- results[grep(pattern = day, results)]
  selected_cluster <- list("TPO" = selected_cluster[grep(pattern = "tpo", selected_cluster)],
                           "PBS" = selected_cluster[grep(pattern = "pbs", selected_cluster)],
                           "fusion" = selected_cluster[grep(pattern = "\\.fusion", selected_cluster)],
                           "nofusion" = selected_cluster[grep(pattern = "nofusion", selected_cluster)])
  
  
  # analyse for each category
  for(category in names(selected_cluster)){
    
    # example RNAseq from HE, differential expression analysis.R lines: position around 2200
    # select the count of genes
    mat_heatmap <- matrix(data = NA, nrow = length(selected_pathways),
                          ncol = length(selected_cluster[[category]]))
    colnames(mat_heatmap) = selected_cluster[[category]]
    rownames(mat_heatmap) = selected_pathways  
    
    
    # update the matrix in the loop
    for(cluster in selected_cluster[[category]]){
      
      # select the temporary table
      tmp_table <- camera_results$h.CATEGORY.v7.5.1.ensembl.mmu[[cluster]]
      
      # remove the unsignificant term
      tmp_table <-  tmp_table[tmp_table$FDR <= 0.05,]
      
      # add for c2 table
      tmp_table <- rbind(tmp_table, camera_results$c2.CATEGORY.v7.5.1.ensembl.mmu[[cluster]])
      
      # remove the unsignificant term
      tmp_table <-  tmp_table[tmp_table$FDR <= 0.05,]
      
      # add information to the table
      if(nrow(tmp_table) > 0){
        
        # select the pathways
        tmp_pathway <- intersect(rownames(tmp_table), rownames(mat_heatmap))
        
        # update the table
        if(length(tmp_pathway) > 0){
          mat_heatmap[tmp_pathway,cluster] = tmp_table[tmp_pathway, ]$Direction  
        }
        
      }
      
    }
    
    
    # change the colnames of the matrix
    colnames(mat_heatmap) = sapply(colnames(mat_heatmap), function(x) strsplit(x, split = "\\.")[[1]][1])
    
    # remove pathways without values
    mat_heatmap = mat_heatmap[!apply(mat_heatmap, 1, function(x) all(is.na(x))), ]
    
    # create the data.frame based on matrix
    dt_heatmap <-  as.data.frame(mat_heatmap)
    dt_heatmap <- stack(dt_heatmap)
    dt_heatmap <- cbind(dt_heatmap, rownames(mat_heatmap))
    
    # change the colnames of the dataframe
    colnames(dt_heatmap) = c("values", "cluster", "pathway")
    
    # add the colours
    dt_heatmap$colour <- "white"
    dt_heatmap$colour[grep(pattern = "Down", dt_heatmap$values)] = "darkgreen"
    dt_heatmap$colour[grep(pattern = "Up", dt_heatmap$values)] = "darkred"
    
    # change the NA colour
    dt_heatmap$values[is.na(dt_heatmap$values)] = "NA"
    
    # update the annotations
    dt_heatmap$cluster = sapply(dt_heatmap$cluster, getAnnotation)
    
    
    # create heatmap with ggplot2
    matplot <- ggplot(dt_heatmap, aes(x = cluster , y = pathway, fill = values)) +
      geom_tile() + theme_bw() +
      theme(axis.text.x = element_text(size = 12, angle = 60, vjust = 0.5),
            axis.text.y = element_text(size = 12),
            axis.title = element_blank(), plot.background = element_rect(fill = "white")) +
      scale_fill_manual(breaks = unique(dt_heatmap$values), values = unique(dt_heatmap$colour))
    
    
    
    # save the figure in the files
    saveFigures(fileName = paste0("matplot_", day, "_", category, "FDR0.05"), ggplot = matplot, 
                dirPlot = outDir, A4 = T)
    

  }
}



# figure 3B -------------------------------------------------------------------------

sce_original <- readRDS("Datasets/sce_day2_and_day5_combined_clutering.rds")

# add the new annotations 
sce_original$final_annotation = as.character(sce_original$combined_clustering)
sce_original$final_annotation[sce_original$combined_clustering %in% c("0.0")] = "ST-HSC_1"
sce_original$final_annotation[sce_original$combined_clustering %in% c("0.2")] = "ST-HSC_2"
sce_original$final_annotation[sce_original$combined_clustering %in% c("0.1")] = "LT-HSC"
sce_original$final_annotation[sce_original$combined_clustering %in% c("1")] = "MPP1_1"
sce_original$final_annotation[sce_original$combined_clustering %in% c("6")] = "MPP1_2"
sce_original$final_annotation[sce_original$combined_clustering %in% c("2")] = "MPP3_1"
sce_original$final_annotation[sce_original$combined_clustering %in% c("4")] = "MPP3_2"
sce_original$final_annotation[sce_original$combined_clustering %in% c("5")] = "MPP3_3"
sce_original$final_annotation[sce_original$combined_clustering %in% c("8")] = "MPP3_4"
sce_original$final_annotation[sce_original$combined_clustering %in% c("3")] = "MPP2_1"
sce_original$final_annotation[sce_original$combined_clustering %in% c("7")] = "MPP2_2"
sce_original$final_annotation[sce_original$combined_clustering %in% c("11")] = "MPP2_3"
sce_original$final_annotation[sce_original$combined_clustering %in% c("12")] = "CMP_1"
sce_original$final_annotation[sce_original$combined_clustering %in% c("16")] = "CMP_2"
sce_original$final_annotation[sce_original$combined_clustering %in% c("9")] = "GMP_1"
sce_original$final_annotation[sce_original$combined_clustering %in% c("10")] = "GMP_2"
sce_original$final_annotation[sce_original$combined_clustering %in% c("13")] = "MEP"
sce_original$final_annotation[sce_original$combined_clustering %in% c("14")] = "T/NK cells pro"
sce_original$final_annotation[sce_original$combined_clustering %in% c("15", "19")] = NA
sce_original$final_annotation[sce_original$combined_clustering %in% c("17")] = "CLP"
sce_original$final_annotation[sce_original$combined_clustering %in% c("18")] = "B cells pro"



# plot the UMAP for the annotations in original sce
plot <- plotReducedDim(sce_original, dimred="TSNE_MNN_Sample_k50", colour_by="final_annotation", point_size = 1) +  fontsize + 
  scale_color_manual(name = "Cell annotation", 
                     breaks = c("B cells pro", "CLP", "T/NK cells pro", "MEP", 
                                "GMP_1", "GMP_2", "CMP_1", "CMP_2", "MPP1_1", 
                                "MPP1_2", "MPP2_1", "MPP2_2", "MPP2_3", "MPP3_1", 
                                "MPP3_2", "MPP3_3", "MPP3_4", "LT-HSC", "ST-HSC_1", "ST-HSC_2"), 
                     values= c('#808000', '#9a6324', '#aaffc3', '#CD00CD',
                               '#0ac5ce', '#03949b', '#7102ff', '#390280', 
                               '#ffe119', '#d4bb15', '#186cf4', '#001afb',
                               '#0214ae', '#3cb44b', '#158022', '#04650f', 
                               '#004808', '#e6194b', '#f58231', '#f85810'))

# save the plot
saveFigures(fileName = "TSNE_cellAnnotation_FinalVersion_fig3", ggplot = plot, dirPlot = outDir)


# figure 3C -------------------------------------------------------------------------


sce <- readRDS("Datasets/sce_iSEE_Outofmemory.rds")

# select the sce cells
sce_original_gating <- sce_original[,colnames(sce)]


# add the gating fast information
sce_original_gating$gating_fast <- sce$gating_fast


# plot the UMAP for the annotations in original sce
plot <- plotReducedDim(sce_original_gating, dimred="TSNE_MNN_Sample_k50", colour_by="gating_fast", point_size = 1) +  fontsize + 
  scale_color_manual(name = "Gating Fast", 
                     breaks = c("LTHSC", "MPP0", "MPP2", "MPP3", 
                                "STHSC", "no_gate"), 
                     values= c('#EB436D', '#E1E435', '#4F64AE', '#40B34B', 
                               '#F4812E', 'gray'))


# save the plot
saveFigures(fileName = "TSNE_gating_FinalVersion_fig3", ggplot = plot, dirPlot = outDir)
