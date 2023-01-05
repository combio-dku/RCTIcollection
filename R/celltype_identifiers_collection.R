## Celltype Identifiers Collection for comparison
## Seokhyun Yoon (syoon@dku.edu), MLBI@DKU, Jan 03, 2023

# suppressPackageStartupMessages(library(igraph))
# suppressPackageStartupMessages(library(scater))
# suppressPackageStartupMessages(library(xgboost))
# suppressPackageStartupMessages(library(SingleCellExperiment))
# suppressPackageStartupMessages(library(dplyr))
# suppressPackageStartupMessages(library(stringr))
# suppressPackageStartupMessages(library(preprocessCore))
# suppressPackageStartupMessages(library(Seurat))
# suppressPackageStartupMessages(library(org.Hs.eg.db))
# suppressPackageStartupMessages(library(scuttle))

# suppressPackageStartupMessages(library(SingleR))
# suppressPackageStartupMessages(library(CHETAH))
# suppressPackageStartupMessages(library(scmap))

# suppressPackageStartupMessages(library(SCINA))
# suppressPackageStartupMessages(library(scSorter))
# suppressPackageStartupMessages(library(garnett))
# suppressPackageStartupMessages(library(scCATCH))

# suppressPackageStartupMessages(library(reticulate))

## Marker-based Identifiers

make_scsorter_ann <- function(df_mkr_mat)
{
  genes <- colnames(df_mkr_mat)
  cell_types <- rownames(df_mkr_mat)

  cnt <- 0
  for(r in 1:dim(df_mkr_mat)[1])
  {
    for(c in 1:dim(df_mkr_mat)[2])
    {
      if(df_mkr_mat[r,c] > 0)
      {
        if(cnt == 0)
        {
          df <- data.frame(cell_types[r], genes[c], 1)
          cnt = cnt + 1
        } else
        {
          df[nrow(df)+1,] <- c(cell_types[r], genes[c], 1)
          cnt = cnt + 1
        }
      }
    }
  }
  colnames(df) <- c('Type','Marker', 'Weight')
  cat(sprintf('scSorter: added %s entries in the marker list\n', cnt))
  flush.console()
  return(df)
}

run_scsorter <- function(df_data, df_mkr_mat)
{
  anno <- make_scsorter_ann(df_mkr_mat)
  expr <- t(df_data)

  cat('scSorter: preprocessing .. ')
  flush.console()
  expr_obj = CreateSeuratObject(counts = expr)
  expr_obj <- NormalizeData(expr_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  expr_obj <- FindVariableFeatures(expr_obj, selection.method = "vst", nfeatures = 2000, verbose = F)
  topgenes <- head(VariableFeatures(expr_obj), 2000)
  expr = GetAssayData(expr_obj)
  # topgene_filter <- rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
  # topgenes <- topgenes[topgene_filter]

  picked_genes = unique(c(anno$Marker, topgenes))
  expr = expr[rownames(expr) %in% picked_genes, ]

  cat('done.\nscSorter: identifying cell types .. ')
  flush.console()
  rts <- scSorter(expr, anno)
  cat('done.\n')
  flush.console()

  return(rts$Pred_Type)
}

# # load libraries and functions
# lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

scType_db_to_mkr_mat <- function(file_name = "cell_markers_ScTypeDB_full.csv",
                                 tissue = c('Immune system'))
{
  df_db <- read.csv(file_name)
  df_db <- df_db[df_db$tissueType == tissue,]

  slst <- list()
  slst2 <- list()
  for(i in 1:nrow(df_db))
  {
    s <- toupper(df_db[i, 'geneSymbolmore1'])
    slst[[i]] <- s
    # print(str_split(s, ',')[[1]])
    if(i == 1){ genes <- str_split(s, ',')[[1]] }
    else{ genes <- union(genes, (str_split(s, ',')[[1]])) }

    s <- toupper(df_db[i, 'geneSymbolmore2'])
    slst2[[i]] <- s
    # print(str_split(s, ',')[[1]])
    if(i == 1){ genes2 <- str_split(s, ',')[[1]] }
    else{ genes2 <- union(genes2, (str_split(s, ',')[[1]])) }

  }
  genes <- union(genes, genes2)
  genes <- sort(genes)
  # print(genes)

  mat <- matrix(rep(0,nrow(df_db)*length(genes)), nrow = nrow(df_db), ncol = length(genes))
  df <- data.frame(mat)
  rownames(df) <- as.vector(df_db[,'cellName'])
  colnames(df) <- genes

  for(i in 1:nrow(df_db))
  {
    genes <- str_split(slst[[i]], ',')[[1]]
    for(g in genes){ df[i,g] <- 1 }

    if(slst2[[i]] != '')
    {
      genes <- str_split(slst2[[i]], ',')[[1]]
      if(length(genes) > 0){for(g in genes){ df[i,g] <- -1 }}
    }
  }

  return(df)
}

mkr_mat_to_scType_df <- function(df_mkr_mat, tissue = 'All')
{
  celltypes <- rownames(df_mkr_mat)
  gene_names <- colnames(df_mkr_mat)
  n_mkrs <- sum(sum(df_mkr_mat > 0))

  df <- data.frame(tissueType = rep(tissue, length(celltypes)))
  df$cellName <- celltypes
  df$geneSymbolmore1 <- rep('', length(celltypes))
  df$geneSymbolmore2 <- rep('', length(celltypes))
  df$shortName <- celltypes

  cnt <- 0
  for(i in 1:nrow(df_mkr_mat))
  {
    b <- df_mkr_mat[i,] > 0
    mkrs <- gene_names[b]
    n <- sum(b)
    s <- ''
    if(length(mkrs) > 0)
    {
      s <- as.character(mkrs[1])
      if(length(mkrs) > 1)
      {
        for(j in 2:length(mkrs))
        {
          s <- sprintf('%s,%s', s, mkrs[j])
        }
      }
    }
    df[i, 'geneSymbolmore1'] <- s
  }
  return(df)
}

scType_df_to_list <- function(df)
{
  cell_types <- df$cellName
  lst_pos <- list()
  lst_neg <- list()
  for(i in 1:length(cell_types))
  {
    v <- str_split(df[i,'geneSymbolmore1'], ',')[[1]]
    if(v[1] != ''){ lst_pos[[i]] <- v }
    else{ lst_pos[[i]] <- list() }

    v <- str_split(df[i,'geneSymbolmore2'], ',')[[1]]
    if(v[1] != ''){ lst_neg[[i]] <- v }
    else{ lst_neg[[i]] <- list() }
  }
  names(lst_pos) <- cell_types
  names(lst_neg) <- cell_types

  return(list(gs_positive = lst_pos, gs_negative = lst_neg))
}

run_sctype <- function(df_data, df_mkr_mat = "cell_markers_ScTypeDB_full.xlsx",
                       Tissues = c('Immune system', 'Pancreas'), N_pca = 15,
                       Clustering_resolution = 1, N_features = 2000)
{
  cat('scType running .. ')
  flush.console()

  if(is.character(df_mkr_mat))
  {
    gs_list = gene_sets_prepare(df_mkr_mat, Tissues)
  } else
  {
    df <- mkr_mat_to_scType_df(df_mkr_mat, tissue = 'All')
    gs_list <- scType_df_to_list(df)
  }

  Seurat_obj <- CreateSeuratObject(t(df_data), project = "SeuratProject", assay = "RNA")
  Seurat_obj[["percent.mt"]] <- PercentageFeatureSet(Seurat_obj, pattern = "^MT-")
  Seurat_obj <- NormalizeData(Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  Seurat_obj <- FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = N_features)
  Seurat_obj <- ScaleData(Seurat_obj, features = rownames(Seurat_obj))
  Seurat_obj <- RunPCA(Seurat_obj, features = VariableFeatures(object = Seurat_obj))
  Seurat_obj <- FindNeighbors(Seurat_obj, dims = 1:N_pca)
  Seurat_obj <- FindClusters(Seurat_obj, resolution = Clustering_resolution)

  es.max = sctype_score(scRNAseqData = Seurat_obj[["RNA"]]@scale.data, scaled = TRUE,
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

  cL_resutls = do.call("rbind", lapply(unique(Seurat_obj@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(Seurat_obj@meta.data[Seurat_obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl,
                    ncells = sum(Seurat_obj@meta.data$seurat_clusters==cl)), 10)}))

  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

  Seurat_obj@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster))
  {
    cl_type = sctype_scores[sctype_scores$cluster==j,];
    Seurat_obj@meta.data$customclassif[Seurat_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }

  cat('scType done.\n')
  flush.console()

  return(Seurat_obj@meta.data$customclassif)
}

scCatch_db_to_mkr_mat <- function(file_name = 'cell_markers_scCatch.csv',
                                  tissue = c('Peripheral blood'), species = 'Human')
{
  if( is.null(file_name) )
  {
    df_db <- cellmatch[cellmatch$species == species & cellmatch$tissue %in% tissue, ]
  } else
  {
    df_db <- read.csv(file_name)
    df_db <- df_db[df_db$species == species & df_db$tissue %in% tissue, ]
  }

  celltypes <- unique(df_db$celltype)
  celltypes <- sort(celltypes)
  genes <- unique(df_db$gene)
  genes <- sort(genes)

  nr <- length(celltypes)
  nc <- length(genes)

  df <- data.frame(matrix(rep(0, nr*nc), nrow = nr, ncol = nc))
  rownames(df) <- celltypes
  colnames(df) <- genes

  for(i in 1:nrow(df_db))
  {
    g <- df_db[i, 'gene']
    c <- df_db[i, 'celltype']
    df[c,g] <- 1
  }

  return(df)
}

mkr_mat_to_scCatch_df <- function(df_mkr_mat, species = 'Human', tissue = 'All')
{
  celltypes <- rownames(df_mkr_mat)
  gene_names <- colnames(df_mkr_mat)
  n_mkrs <- sum(sum(df_mkr_mat > 0))

  df <- data.frame(species = rep(species, n_mkrs), row.names = seq(n_mkrs))
  df$tissue <- rep(tissue, n_mkrs)
  df$cancer <- rep('Normal', n_mkrs)
  df$condition <- rep('Normal cell', n_mkrs)
  df$subtype1 <- rep('NA', n_mkrs)
  df$subtype2 <- rep('NA', n_mkrs)
  df$subtype3 <- rep('NA', n_mkrs)
  df$celltype <- rep('NA', n_mkrs)
  df$gene <- rep('NA', n_mkrs)
  df$resource <- rep('Experiment', n_mkrs)
  df$pmid <- rep('NA', n_mkrs)

  cnt <- 0
  for(i in 1:nrow(df_mkr_mat))
  {
    b <- df_mkr_mat[i,] > 0
    mkrs <- gene_names[b]
    n <- sum(b)
    df[(cnt+1):(cnt+n), 'celltype'] <- celltypes[i]
    df[(cnt+1):(cnt+n), 'gene'] <- mkrs
  }
  return(df)
}

run_sccatch <- function(df_data, df_mkr_mat = 'cell_markers_scCatch.csv',
                        Tissues = c(), Species = 'Human',
                        N_pca = 15, Clustering_resolution = 1, N_features = 2000)
{
  cat('scCatch running .. ')
  flush.console()

  Seurat_obj<-CreateSeuratObject(t(df_data), assay = "RNA")
  Seurat_obj[["percent.mt"]] <- PercentageFeatureSet(Seurat_obj, pattern = "^MT-")
  Seurat_obj <- NormalizeData(Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  Seurat_obj <- FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = N_features)
  Seurat_obj <- ScaleData(Seurat_obj, features = rownames(Seurat_obj))
  Seurat_obj <- RunPCA(Seurat_obj, features = VariableFeatures(object = Seurat_obj))
  Seurat_obj <- FindNeighbors(Seurat_obj, dims = 1:N_pca)
  Seurat_obj <- FindClusters(Seurat_obj, resolution = Clustering_resolution)

  Seurat_obj@meta.data$clusters <- as.character(Seurat_obj@meta.data$seurat_clusters)
  catch_obj <- createscCATCH(Seurat_obj@assays$RNA@data,
                             cluster = Seurat_obj@meta.data$clusters)

  if( is.null(df_mkr_mat) )
  {
    cellmatch_new <- cellmatch[cellmatch$species == Species & cellmatch$tissue %in% Tissues, ]
  } else if( is.character(df_mkr_mat) )
  {
    cellmatch_new <- read.csv(df_mkr_mat)
    cellmatch_new <- cellmatch_new[cellmatch_new$species == Species & cellmatch_new$tissue %in% Tissues, ]
  } else
  {
    cellmatch_new <- mkr_mat_to_scCatch_df(df_mkr_mat, species = Species, tissue = 'All')
  }
  head(cellmatch_new)

  catch_obj <- findmarkergene(object = catch_obj, if_use_custom_marker = TRUE,
                              marker = cellmatch_new)
  catch_obj <- findcelltype(object = catch_obj)

  Seurat_obj@meta.data[,'scCatch_pred'] <- 'Unknown'
  for(j in 1:nrow(catch_obj@celltype))
  {
    b <- Seurat_obj@meta.data$clusters == catch_obj@celltype[j,'cluster']
    Seurat_obj@meta.data[b,'scCatch_pred'] <- catch_obj@celltype[j,'cell_type']
  }

  cat('scCatch done.\n')
  flush.console()

  return(Seurat_obj@meta.data$scCatch_pred)
}

make_scina_ann <- function(df_mkr_mat, genes_data)
{
  genes <- colnames(df_mkr_mat)
  cell_types <- rownames(df_mkr_mat)

  lst <- list()
  for(r in 1:dim(df_mkr_mat)[1])
  {
    genes_sel <- genes[df_mkr_mat[r,] > 0]
    lst <- c(lst, list(genes_sel))
    x <- intersect(genes_sel, genes_data)
    if(length(x) == 0)
    {
      cat(sprintf('No genes for %s', cell_types[r]))
    }
  }
  names(lst) <- cell_types

  return(lst)
}

run_scina <- function(df_data, df_mkr_mat, s_cutoff = 0.9)
{
  cat('SCINA running .. ')
  flush.console()
  signatures <- make_scina_ann(df_mkr_mat, colnames(df_data))

  exp_raw <- t(df_data)
  exp_raw <- log(exp_raw+1)
  exp_raw[] <- (normalize.quantiles(exp_raw))

  results = SCINA(exp_raw, signatures, max_iter = 100, convergence_n = 10,
                  convergence_rate = 0.999, sensitivity_cutoff = s_cutoff,
                  rm_overlap=FALSE, allow_unknown=TRUE, log_file='SCINA.log')
  cat('done.\n')
  flush.console()

  return(results$cell_labels)
}


gen_garnett_marker_file_old <- function(df_mkr_mat, file_name = 'Garnett_markers.txt')
{
  sink(file_name)
  ct_lst <- rownames(df_mkr_mat)
  genes <- colnames(df_mkr_mat)

  for(i in 1:length(ct_lst))
  {
    ct <- ct_lst[i]
    b <- (df_mkr_mat[ct,] > 0)
    gsel <- genes[b]

    if( sum(b) > 0)
    {
      cat(paste0('>',ct,'\nexpressed: '))
      for(m in 1:length(gsel))
      {
        cat(gsel[m])
        if(m < length(gsel)){cat(', ')}
      }
      cat('\n')
    }

  }
  sink()
}

gen_garnett_marker_file <- function(df_mkr_mat, file_name = 'Garnett_markers.txt')
{

  ct_lst <- rownames(df_mkr_mat)
  genes <- colnames(df_mkr_mat)

  lines <- ''
  for(i in 1:length(ct_lst))
  {
    ct <- ct_lst[i]
    b <- (df_mkr_mat[ct,] > 0)
    gsel <- genes[b]

    if( sum(b) > 0)
    {
      s <- paste0('>',ct,'\nexpressed: ')
      for(m in 1:length(gsel))
      {
        s <- paste0( s, sprintf('%s', gsel[m]) )
        if(m < length(gsel))
        {
          s <- paste0( s, ', ' )
        }
      }
      s <- paste0( s, '\n' )
    }
    lines <- paste0(lines, s)
  }
  writeLines(lines, file_name)
}

pbmc_cell_type_map <- list( 'NK cells' = 'NK cell',
                            'Monocytes' = 'Monocyte',
                            'B cells' = 'B cell',
                            'T cells' = 'T cell',
                            'CD4 T cells' = 'T cell',
                            'CD8 T cells' = 'T cell',
                            'Dendritic cells' = 'Dendritic cell')

rename_cell_type_garnet <- function( cell_type_org, ct_map_lst)
{
  ct_lst <- names(ct_map_lst)
  for(i in 1:length(ct_lst))
  {
    cell_type_org[cell_type_org == ct_lst[i]] <- ct_map_lst[[ct_lst[i]]]

  }
  return(cell_type_org)
}

## Works only for PBMC
run_garnett <- function(df_data, df_mkr_mat, classifier = NA)
{
  cat('Run Garnett .. ')
  flush.console()
  mkr_file <- 'Garnett_markers_tmp.txt'
  gen_garnett_marker_file(df_mkr_mat, file_name = mkr_file)

  # load in the data
  # NOTE: the 'system.file' file name is only necessary to read in
  # included package data
  #
  mat <- as( t(df_data), "dgCMatrix")
  pdata <- data.frame( cell = colnames(mat) )
  fdata <- data.frame( gene_short_name = rownames(mat), row.names = rownames(mat) )

  row.names(mat) <- row.names(fdata)
  colnames(mat) <- row.names(pdata)

  # create a new CDS object
  pd <- new("AnnotatedDataFrame", data = pdata)
  fd <- new("AnnotatedDataFrame", data = fdata)
  pbmc_cds <- newCellDataSet(as(mat, "dgCMatrix"),
                             phenoData = pd,
                             featureData = fd)

  # generate size factors for normalization later
  pbmc_cds <- estimateSizeFactors(pbmc_cds)

  classifier <- train_cell_classifier(cds = pbmc_cds,
                                      marker_file = mkr_file,
                                      db=org.Hs.eg.db,
                                      cds_gene_id_type = "SYMBOL",
                                      num_unknown = 50,
                                      marker_file_gene_id_type = "SYMBOL")

  # if(is.na(classifier))
  # {
  #     classifier <- readRDS("hsPBMC_20191017.RDS")
  # }

  pbmc_cds <- classify_cells(pbmc_cds, classifier,
                             db = org.Hs.eg.db,
                             cluster_extend = TRUE,
                             cds_gene_id_type = "SYMBOL")


  # head(pData(pbmc_cds))
  # table(pData(pbmc_cds)$cell_type)
  # table(pData(pbmc_cds)$cluster_ext_type)

  pdata[,'cell_type_pred'] <- pData(pbmc_cds)[,'cell_type']
  # pdata[,'cell_type_pred'] <- rename_cell_type_garnet( pdata[,'cell_type_pred'], pbmc_cell_type_map)

  cat('done.\n')
  flush.console()

  return(pdata[,'cell_type_pred'])
}


run_hicat <- function( df_data, marker_file = 'cell_markers_rndsystems_rev.tsv',
                       Tissues = c(), N_pca = 15, Clustering_resolution = 1 )
{
  if( !is.character(marker_file) )
  {
    cat('ERROR: Marker must be a string representing the full path to the marker file to use.')
    return(NULL)
  }
  ### Load MarkerCount
  hicat <- import("MarkerCount.hicat")

  cat('Run HiCAT .. ')
  flush.console()

  lst_res <- hicat$HiCAT( df_data, marker_file = marker_file, log_transformed = FALSE,
                          N_pca_components = N_pca, # target_tissues = Tissues,
                          Clustering_resolution = Clustering_resolution )

  df_pred <- lst_res[[1]]
  summary <- lst_res[[2]]

  cat('done.\n')
  flush.console()

  return(df_pred)
}

#' identify_celltypes_using_marker
#'
#' This function runs one of the marker-based celltype identifiers, specified by method.
#'
#' @param df_data input cell-by-gene matrix (data frame), column name must be in Hugo symbol
#' @param marker it can be either a string or a data frame (marker-matrix). For 'HiCAT', it must be a string representing the full path to the HiCAT marker file. For others, it must be a data frame obtained by running 'get_mkr_mat_from_hicat_mkr_file'. 
#' @param method method to use. It can be one of 'HiCAT', 'SCINA', 'scSorter', 'scType', 'scCatch', or 'Garnett'.
#' @return a vector of predicted celltypes. Each element corresponds to a row of df_data. For HiCAT, it is a data frame with 3 columns, each corresponds to major-type, minor-type and subset.
#' @export
identify_celltypes_using_marker <- function( df_data, marker, method = NULL,
                                             N_pca = 15,  Clustering_resolution = 1, N_features = 2000)
{
  df_mkr <- marker
  if( is.character(df_mkr) )
  {
    ## pass
  } else
  {
    genes_d <- colnames(df_data)
    genes_c <- colnames(df_mkr)

    genes <- intersect(genes_d, genes_c)
    df_mkr <- df_mkr[,genes]

    ctypes_lst <- rownames(df_mkr)
  }

  if(tolower(method) == 'scina')
  {
    celltypes <- run_scina(df_data, df_mkr, s_cutoff = 0.9)
  } else if(tolower(method) == 'scsorter')
  {
    celltypes <- run_scsorter(df_data, df_mkr)
  } else if(tolower(method) == 'garnett')
  {
    celltypes <- run_garnett(df_data, df_mkr, classifier = NA)
  } else if(tolower(method) == 'sccatch')
  {
    celltypes <- run_sccatch(df_data, df_mkr_mat = df_mkr,
                             Tissues = c(), Species = 'Human',
                             N_pca = N_pca,
                             Clustering_resolution = Clustering_resolution,
                             N_features = N_features)
  } else if(tolower(method) == 'sctype')
  {
    celltypes <- run_sctype(df_data, df_mkr_mat = df_mkr,
                            Tissues = c(),
                            N_pca = N_pca,
                            Clustering_resolution = Clustering_resolution,
                            N_features = N_features)
  } else if(tolower(method) == 'hicat')
  {
    celltypes <- run_hicat( df_data, marker_file = df_mkr,
                            Tissues = c(),
                            N_pca = N_pca,
                            Clustering_resolution = Clustering_resolution )
  } else
  {
    s <- sprintf('method not specified. \npossible method: scina, sctype, scsorter, garnett or sccatch\n')
    cat(s)
  }

  return(celltypes)
}


## Reference-based Identifiers

run_singler <- function( df_data_test, df_data_ref, celltypes_ref,
                         log_transformed = FALSE, undesired_lst = list() )
{
  cat('Run SingleR .. ')
  flush.console()
  df_ct_ref <- data.frame( cell_type = celltypes_ref )
  rownames(df_ct_ref) <- rownames(df_data_ref)

  ## Make SingleCellExperiments
  reference <- SingleCellExperiment(assays = list(counts = t(df_data_ref)),
                                    colData = df_ct_ref )

  input <- SingleCellExperiment(assays = list(counts = t(df_data_test)))

  if(!log_transformed)
  {
    reference <- logNormCounts(reference)
    input <- logNormCounts(input)
  }

  ## Run CHETAH
  pred.grun <- SingleR(test=input, ref=reference, labels=df_ct_ref[,'cell_type'], de.method="wilcox")

  cat('done. \n')
  ## Extract celltypes:
  # celltypes <- pred.grun$labels

  return(pred.grun$labels)
}

run_chetah <- function( df_data_test, df_data_ref, celltypes_ref,
                        log_transformed = FALSE, undesired_lst = list() )
{
  cat('Run CHETAH .. ')
  flush.console()

  df_ct_ref <- data.frame( celltypes = celltypes_ref )
  rownames(df_ct_ref) <- rownames(df_data_ref)

  if(!log_transformed)
  {
    df_data_ref <- log2(1+df_data_ref)
    df_data_test <- log2(1+df_data_test)
  }

  ## Make SingleCellExperiments
  reference <- SingleCellExperiment(assays = list(counts = t(df_data_ref)),
                                    colData = df_ct_ref )

  input <- SingleCellExperiment(assays = list(counts = t(df_data_test)))

  if(!log_transformed)
  {
    reference <- logNormCounts(reference)
    # logcounts(reference) <- log2(normcounts(reference) + 1)

    input <- logNormCounts(input)
    # logcounts(input) <- log2(normcounts(input) + 1)
  }


  ## Run CHETAH
  input <- CHETAHclassifier(input = input, ref_cells = reference)

  ## Extract celltypes:
  # celltypes <- input$celltype_CHETAH

  cat('done.\n')
  flush.console()

  return(input$celltype_CHETAH)
}

run_scmap <- function( df_data_test, df_data_ref, celltypes_ref,
                       log_transformed = FALSE, undesired_lst = list(), sel = 'cell' )
{
  ref_mat <- t(df_data_ref) # row: genes, column: cell
  cell_type <- data.frame( cell_type1 = celltypes_ref )
  rownames(cell_type) <- rownames(df_data_ref)

  cat(sprintf('scmap_%s: preprocessing data .. ', sel))
  flush.console()
  ref <- SingleCellExperiment(assays = list(counts = as.matrix(ref_mat)), colData = cell_type)

  test_mat <- t(df_data_test)
  test <- SingleCellExperiment(assays = list(counts = as.matrix(test_mat)))

  if(!log_transformed)
  {
    ref <- logNormCounts(ref)
    test <- logNormCounts(test)
    # logcounts(ref) <- log2(normcounts(ref) + 1)
    # logcounts(test) <- log2(normcounts(test) + 1)
  } else
  {
    # logcounts(ref) <- normcounts(ref)
  }
  # use gene names as feature symbols
  rowData(ref)$feature_symbol <- rownames(ref)
  # remove features with duplicated names
  ref <- ref[!duplicated(rownames(ref)), ]

  # use gene names as feature symbols
  rowData(test)$feature_symbol <- rownames(test)
  # remove features with duplicated names
  test <- test[!duplicated(rownames(test)), ]

  cat(sprintf('done\nscmap_%s: preparing reference .. ', sel))
  flush.console()
  ## Select features
  ref <- selectFeatures(ref, suppress_plot = TRUE)

  if(sel == 'cluster')
  {
    ref <- indexCluster(ref)

    cat(sprintf('done\nscmap_%s: identifying cell types .. ', sel))
    flush.console()
    scmap_results <- scmapCluster(projection = test,
                                  index_list = list(
                                    cell_type_pred = metadata(ref)$scmap_cluster_index
                                  )
    )
  } else
  {
    ref <- indexCell(ref)

    cat(sprintf('done\nscmap_%s: identifying cell types .. ', sel))
    flush.console()
    scmapCell_results <- scmapCell(
      test,
      list(
        cell_type_pred = metadata(ref)$scmap_cell_index
      )
    )

    scmap_results <- scmapCell2Cluster(
      scmapCell_results,
      list(
        as.character(colData(ref)$cell_type1)
      )
    )
  }
  cat('done\n')

  return( scmap_results$combined_labs)
}

run_castle <- function(df_data_test, df_data_ref, celltypes_ref,
                       log_transformed = FALSE, undesired_lst = list(),
                       UA_threshold = 0.5, nFeatures = 1000, BREAKS=c(-1, 0, 1, 6, Inf))
{
  if( log_transformed )
  {
    cat('CaSTLe requires raw count matrix, not log-normalized GEP.\n')
    cat('Use raw count matrix with "log_transform = FALSE" \n')
    return(NULL)
  }

  df_ct_ref <- data.frame( celltypes = celltypes_ref )
  rownames(df_ct_ref) <- rownames(df_data_ref)

  ds1 <- df_data_ref
  ds2 <- df_data_test
  sourceCellTypes <- as.factor(df_ct_ref[,'celltypes'])

  cat('CaSTLe preprocessing .. ')
  flush.console()
  source <- SingleCellExperiment(assays = list(counts = t(df_data_ref)),
                                 colData = df_ct_ref )

  target <- SingleCellExperiment(assays = list(counts = t(df_data_test)))

  # remove unlabeled from source
  unknownLabels = levels(sourceCellTypes)[grep("not applicable|unclassified|contaminated|unknown",
                                               levels(sourceCellTypes))]
  if (length(unknownLabels)>0) {
    hasKnownLabel = is.na(match(sourceCellTypes, unknownLabels))
    sourceCellTypes = sourceCellTypes[hasKnownLabel]
    sourceCellTypes = as.factor(as.character(sourceCellTypes))
    ds1 = ds1[hasKnownLabel,]
  }

  # 2. Unify sets, excluding low expressed genes
  source_n_cells_counts = apply(counts(source), 1, function(x) { sum(x > 0) } )
  target_n_cells_counts = apply(counts(target), 1, function(x) { sum(x > 0) } )
  common_genes = intersect( colnames(ds1)[source_n_cells_counts>10],
                            colnames(ds2)[target_n_cells_counts>10]
  )
  remove(source_n_cells_counts, target_n_cells_counts)
  ds1 = ds1[, colnames(ds1) %in% common_genes]
  ds2 = ds2[, colnames(ds2) %in% common_genes]
  ds = rbind(ds1[,common_genes], ds2[,common_genes])
  isSource = c(rep(TRUE,nrow(ds1)), rep(FALSE,nrow(ds2)))
  remove(ds1, ds2)

  # 3. Highest mean in both source and target
  topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]

  # for each cell - what is the most probable classification?
  L = length(levels(sourceCellTypes))
  targetClassification = as.data.frame(matrix(rep(0,L*sum(!isSource)), nrow=L),
                                       row.names = levels(sourceCellTypes))

  cat('done.\nCaSTLe cell profiling .. ')
  flush.console()

  # iterate over all source cell types
  for (cellType in levels(sourceCellTypes)) {

    inSourceCellType = as.factor(ifelse(sourceCellTypes == cellType, cellType, paste0("NOT",cellType)))

    # 4. Highest mutual information in source
    topFeaturesMi = names(sort(apply(ds[isSource,],2,
                                     function(x) { compare(cut(x,breaks=BREAKS),inSourceCellType,method = "nmi") }),
                               decreasing = T))

    # 5. Top n genes that appear in both mi and avg
    selectedFeatures = union(head(topFeaturesAvg, nFeatures) , head(topFeaturesMi, nFeatures) )

    # 6. remove correlated features
    tmp = cor(ds[,selectedFeatures], method = "pearson")
    tmp[!lower.tri(tmp)] = 0
    selectedFeatures = selectedFeatures[apply(tmp,2,function(x) any(x < 0.9))]
    remove(tmp)

    # 7,8. Convert data from continous to binned dummy vars
    # break datasets to bins
    dsBins = apply(ds[, selectedFeatures], 2, cut, breaks= BREAKS)
    # use only bins with more than one value
    nUniq = apply(dsBins, 2, function(x) { length(unique(x)) })
    # convert to dummy vars
    ds0 = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))
    remove(dsBins, nUniq)

    cat(paste0(cellType,", "))
    flush.console()

    inTypeSource = sourceCellTypes == cellType
    # 9. Classify
    xg=xgboost(data=ds0[isSource,] ,
               label=inTypeSource,
               objective="binary:logistic",
               eta=0.7 , nthread=1, nround=20, verbose=0,
               gamma=0.001, max_depth=5, min_child_weight=10,
               eval_metric = 'logloss')

    # 10. Predict
    inTypeProb = predict(xg, ds0[!isSource, ])

    targetClassification[cellType,] = inTypeProb
  }
  cat('done.\n')
  flush.console()

  df <- data.frame( t(targetClassification) )
  df <- df/rowSums(df)
  pred_label <- (names(df)[-1][max.col(df[-1], 'first')])
  max_val <- rowMaxs(as.matrix(df))
  pred_label[max_val < UA_threshold] <- 'unassigned'

  for(k in 1:length(pred_label))
  {
    pred_label[k] <- str_replace(pred_label[k], c('[.]'), c(' '))
  }

  return(pred_label)
}


run_markercount_ref <- function( df_data_test, df_data_ref, celltypes_ref,
                                 log_transformed = FALSE, undesired_lst = list() )
{
  ### Load MarkerCount
  mkrcnt <- import("MarkerCount.marker_count")

  cat('Run MarkerCount_Ref .. ')
  flush.console()

  df_ct_ref <- data.frame( celltypes = celltypes_ref )
  rownames(df_ct_ref) <- rownames(df_data_ref)

  df_res <- mkrcnt$MarkerCount_Ref( df_data_ref, celltypes_ref,
                                    X_test = df_data_test, df_mkr_mat = NULL,
                                    N_mkrs = 18, of_th = 0.8, # min_th = 0.2,
                                    cell_types_to_excl = undesired_lst,
                                    cluster_label = NULL, X_pca = NULL,
                                    log_transformed = log_transformed,
                                    # file_to_save_marker = 'UC_cell_markers',
                                    verbose = TRUE)

  cat('done.\n')
  flush.console()

  return(df_res$cell_type_pred)
}


remove_undesired_cell_type <- function(dfd, celltypes, undesired_lst)
{
  bu = rep(FALSE, length(celltypes))
  for(i in 1:length(undesired_lst))
  {
    bu <- bu | (celltypes == undesired_lst[i])
  }
  dfd <- data.frame(dfd[!bu,])
  celltypes <- celltypes[!bu]

  out <- list(dfd, celltypes)
  return(out)
}

#' identify_celltypes_using_ref
#'
#' This function runs one of the reference-based celltype identifiers, specified by method.
#'
#' @param df_data_test input cell-by-gene matrix (data frame) to test. column name (marker gene) must be in Hugo symbol.
#' @param df_data_ref cell-by-gene matrix (data frame) to be used as reference
#' @param celltypes_ref a vector of celltypes corresponding to each row of df_data_ref
#' @param method method to use. It can be one of 'MarkerCount', 'SingleR', 'CHETAH', 'CaSTLe', 'scmap_cell', or 'scmap_cluster'.
#' @param log_transformed If df_data_test is log-normalized, you should set it TRUE. For now, it must be FALSE.
#' @param undesired_lst a list (vector) of celltypes that you don't want to include in reference celltype, e.g., 'Unknown' or 'Tumor'.
#' @return a vector of predicted celltypes. Each element corresponds to a row of df_data_test.
#' @export
identify_celltypes_using_ref <- function( df_data_test, df_data_ref, celltypes_ref, method = NULL,
                                          log_transformed = FALSE, undesired_lst = c() )
{
  if(length(undesired_lst) > 0)
  {
    out <- remove_undesired_cell_type(df_data_ref, celltypes_ref, undesired_lst)
    df_data_ref <- out[[1]]
    celltypes_ref <-  out[[2]]
  }

  if(tolower(method) == 'castle')
  {
    celltypes <- run_castle( df_data_test, df_data_ref, celltypes_ref,
                             log_transformed = log_transformed, UA_threshold = 0.5, undesired_lst )
  } else if(tolower(method) == 'singler')
  {
    celltypes <- run_singler( df_data_test, df_data_ref, celltypes_ref,
                              log_transformed = log_transformed, undesired_lst )
  } else if(tolower(method) == 'scmap_cell')
  {
    celltypes <- run_scmap( df_data_test, df_data_ref, celltypes_ref,
                            log_transformed = log_transformed, sel = 'cell', undesired_lst )
  } else if(tolower(method) == 'scmap_cluster')
  {
    celltypes <- run_scmap( df_data_test, df_data_ref, celltypes_ref,
                            log_transformed = log_transformed, sel = 'cluster', undesired_lst )
  } else if(tolower(method) == 'chetah')
  {
    celltypes <- run_chetah( df_data_test, df_data_ref, celltypes_ref,
                             log_transformed = log_transformed, undesired_lst )
  } else if((tolower(method) == 'markercount') | (tolower(method) == 'markercount_ref'))
  {
    celltypes <- run_markercount_ref( df_data_test, df_data_ref, celltypes_ref,
                                      log_transformed = log_transformed, undesired_lst )
  } else
  {
    s <- sprintf('method not specified. \n')
    cat(s)
    s <- sprintf('supported method: singler, chetah, castle, scma_cell or scmap_cluster\n')
    cat(s)
  }

  return(celltypes)
}

#' Get marker-matrix from HiCAT marker file
#'
#' This function uses HiCAT marker file to generate marker-matrix to be used for
#' celltype identifiers other than HiCAT. Marker-matrix is a data frame with its
#' rows being celltypes and columns being markers. Each entry is 0 or 1 indicating
#' the inclusion of the marker in a celltype.
#'
#' @param mkr_file full path to the HiCAT marker file.
#' @param taxo_level taxonomy level. It can be one of 'major', 'minor' or 'subset'.
#' @param tissues a vector of tissues (groups of celltypes) to be exported. Use tissue names in the HiCAT marker file. 
#' @return marker-matrix in data frame
#' @export
get_mkr_mat_from_hicat_mkr_file <- function( mkr_file, taxo_level = 'major', tissues = NULL)
{
  df <- read.csv(mkr_file, header = TRUE, sep = '\t')
  rns <- rownames(df)

  ct_col <- paste0('cell_type_', taxo_level)

  if(!is.null(tissues))
  {
    b <- df[,'tissue'] %in% tissues
    df <- df[b, ]
    rns <- rownames(df)
  }

  ct_lst <- unique(df[,ct_col])
  mkrs_lst <- list()
  mkrs_all <- c()

  for(k in 1:length(ct_lst))
  {
    ct <- ct_lst[k]

    b1 <- df[,ct_col] == ct
    b2 <- df[,'exp'] == 'pos'
    b <- b1 & b2
    rns_sel <- rns[b]

    mkrs <- c()
    cnt <- 0
    for(rn in rns_sel)
    {
      gs <- df[rn, 'markers']
      genes <- str_split(gs, ',')[[1]]
      # genes_new <- c()
      for(i in 1:length(genes))
      {
        cnt <- cnt + 1
        mkrs[cnt] <- str_trim(genes[i], side = 'both')
      }
    }
    mkrs <- unique(mkrs)
    mkrs_lst[[k]] <- mkrs

    mkrs_all <- union(mkrs_all, mkrs)
    # print(mkrs)
  }

  names(mkrs_lst) <- ct_lst

  mkrs_all <- unique(mkrs_all)
  mkrs_all <- sort(mkrs_all)

  mkr_mat <- matrix(0, length(ct_lst), length(mkrs_all))
  df_mkr_mat <- data.frame(mkr_mat)
  rownames(df_mkr_mat) <- ct_lst
  colnames(df_mkr_mat) <- mkrs_all

  for(k in 1:length(ct_lst))
  {
    ct <- ct_lst[k]
    mkrs <- mkrs_lst[[k]]
    df_mkr_mat[ct, mkrs] <- 1
  }
  return(df_mkr_mat)
}
