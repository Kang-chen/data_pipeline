library(tidyverse)
library(scfetch)
library(Seurat)


library(logger)
log_threshold(TRACE)
log_formatter(formatter_glue_or_sprintf)
# log_layout(layout_glue_colors)

source("/home/rstudio/data_pipeline/R/scfetch/R/GEO.R")
source("/home/rstudio/data_pipeline/R/scfetch/R/utils.R")

Sys.setenv("VROOM_CONNECTION_SIZE"=131072*60)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Please provide the path to the JSON file as a command line argument.")
}

source_id <- args[1]
df_dataset = readxl::read_xlsx("Data_collection.xlsx")

df_query = df_dataset %>% filter(source_id == !!source_id)
log_info(df_query)
file_type = df_query$file_type # rdata
meta_file_pattern = df_query$meta_file_pattern
manual.download =  df_query$manual_download_dir
if (is.na(manual.download)) {
  manual.download = NULL
}

workspace = file.path(getwd(),source_id)
if (!dir.exists(workspace)) {
  dir.create(workspace, recursive = TRUE)
  log_info("Directory created: ", workspace)
}



seurat <- ParseGEO(
  acce = source_id, platform = NULL, supp.idx = 1, down.supp = TRUE, supp.type = file_type,
  out.folder =workspace, gene2feature = F,manual.download=manual.download
)

seurat %>% checkSeuratAssay()

meta <- tryCatch(
  {
    ExtractGEOMetaFile(acce = source_id, 
                       workspace = workspace,
                       meta_file_pattern = meta_file_pattern)
  },
  error = function(e) {
    log_warn("Error in extracting meta data.")
    return(NULL)  
  }
)

if (!is.null(rownames(meta))) {
  if (all(rownames(meta) %in% colnames(seurat))) {
    seurat <- AddMetaData(seurat, metadata = meta)
    log_info("Meta data has been added to the Seurat object.")
  } else {
    warning("The rownames of meta do not match the colnames of the Seurat object.")
    rownames(meta) %>% 
      setdiff(colnames(seurat)) %>% 
      head() %>% 
      log_warn()
  }
}else{
  log_warn("no meta data import!")
}



seurat %>% ExportSeurat(to = "AnnData",conda.path = "/opt/conda",anndata.file = file.path(".",source_id,"outs/processed.h5ad"))