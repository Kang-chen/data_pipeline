#' Download Matrix from GEO and Load to Seurat/DESeq2.
#'
#' @param acce GEO accession number.
#' @param platform Platform information/field. Disable when \code{down.supp} is TRUE. Default: NULL (disable).
#' @param down.supp Logical value, whether to download supplementary files to create count matrix. If TRUE, always
#' download supplementary files. If FALSE, use \code{ExpressionSet} (If contains non-integer or empty,
#' download supplementary files automatically). Default: FALSE.
#' @param supp.idx The index of supplementary files to download. This should be consistent with \code{platform}. Default: 1.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param data.type The data type of the dataset, choose from "sc" (single-cell) and "bulk" (bulk). Default: "sc".
#' @param supp.type The type of downloaded supplementary files, choose from count (count matrix file or single count matrix file),
#' 10x (cellranger output files in tar/gz supplementary files, contains barcodes, genes/features and matrix, e.g. GSE200257)
#' and 10xSingle (cellranger output files in supplementary files directly, e.g. GSE236082). Default: count.
#' @param out.folder Output folder to save 10x files. Default: NULL (current working directory).
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}. Default: TRUE.
#' @param merge Logical value, whether to merge Seurat list when there are multiple 10x files (\code{supp.type} is 10x). Default: FALSE.
#' @param meta.data Dataframe contains sample information for DESeqDataSet, use when \code{data.type} is bulk. Default: NULL.
#' @param fmu Column of \code{meta.data} contains group information. Default: NULL.
#' @param ... Parameters for \code{\link{getGEO}}.
#'
#' @return If \code{data.type} is "sc", return Seurat object (if \code{merge} is TRUE) or Seurat object list (if \code{merge} is FALSE).
#' If \code{data.type} is "bulk", return DESeqDataSet.
#' @importFrom magrittr %>%
#' @importFrom GEOquery getGEO getGEOSuppFiles gunzip
#' @importFrom Biobase annotation experimentData pData phenoData notes sampleNames exprs
#' @importFrom tools file_ext
#' @importFrom utils untar
#' @importFrom data.table fread
#' @importFrom openxlsx read.xlsx
#' @importFrom Seurat Read10X CreateSeuratObject
#' @importFrom methods new
#' @importFrom stats formula
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @export
#'
#' @examples
#' \dontrun{
#' # the supp files are count matrix
#' GSE94820.seu <- ParseGEO(acce = "GSE94820", down.supp = TRUE, supp.idx = 1, supp.type = "count")
#' # the supp files are cellranger output files: barcodes, genes/features and matrix
#' # need users to provide the output folder
#' GSE200257.seu <- ParseGEO(
#'   acce = "GSE200257", down.supp = TRUE, supp.idx = 1, supp.type = "10x",
#'   out.folder = "/path/to/output/folder"
#' )
#' # need users to provide the output folder
#' GSE236082.seu <- ParseGEO(
#'   acce = "GSE236082", down.supp = TRUE, supp.type = "10xSingle",
#'   out.folder = "/path/to/output/folder"
#' )
#' }
ParseGEO <- function(acce, platform = NULL, down.supp = FALSE, supp.idx = 1, timeout = 3600,
                     data.type = c("sc", "bulk"),
                     supp.type = c("count", "10x", "10xSingle","Seurat","SingleCellExperiment","Loom","AnnData"),
                     out.folder = NULL, gene2feature = TRUE, merge = TRUE,
                     meta.data = NULL, fmu = NULL,manual.download = NULL, ...) {
  # check parameters
  data.type <- match.arg(arg = data.type)
  supp.type <- match.arg(arg = supp.type)
  
  if (supp.type %in% c("Seurat","SingleCellExperiment","Loom","AnnData")) {
    ExtractGEOExpSuppObj(acce, timeout = 3600000, supp.idx = 1,manual.download = manual.download,
                                     supp.type=c("Seurat","SingleCellExperiment","Loom","AnnData")) 
  }else{
  
    # check platform
    if (down.supp) {
      message("Download supplementary files to generate matrix!")
      pf.obj <- NULL
    } else {
      message("Extract expression data from eSets!")
      if (is.null(platform)) {
        stop("Platform is required to extract expression data!")
      }
      # get GEO object
      pf.obj <- GEOobj(acce = acce, platform = platform, ...)
    }
    # change supp type to count when bulk
    if (data.type == "bulk") {
      supp.type <- "count"
    }
    # extract counts matrix
    pf.count <- ExtractGEOExp(
      pf.obj = pf.obj, acce = acce, supp.idx = supp.idx, down.supp = down.supp,manual.download= manual.download,
      timeout = timeout, supp.type = supp.type, out.folder = out.folder, gene2feature = gene2feature
    )
    if (data.type == "bulk") {
      de.obj <- Loading2DESeq2(mat = pf.count, meta = meta.data, fmu = fmu)
      return(de.obj)
    } else if (data.type == "sc") {
      # load seurat
      if (is.null(pf.count) && (supp.type == "10x" || supp.type == "10xSingle")) {
        message("Loading data to Seurat!")
        all.samples.folder <- dir(out.folder, full.names = TRUE)
        # check file
        valid.samples.folder <- Check10XFiles(folders = all.samples.folder, gene2feature = gene2feature)
        if (length(valid.samples.folder) == 0) {
          stop("No valid sample folder detected under ", out.folder, ". Please check!")
        }
        # load to seurat
        seu.list <- sapply(valid.samples.folder, function(x) {
          seu.obj <- tryCatch({
            # First try using Read10X_GEO
            if(file.exists(file.path(x,"barcodes.tsv.gz"))){
              checkBarcodes(file.path(x,"barcodes.tsv.gz"))
            }
            x.mat <- scCustomize::Read10X_GEO(data_dir = x)
            Seurat::CreateSeuratObject(counts = x.mat[[1]], project = basename(x))
          }, error = function(e) {
            tryCatch({
              # If Read10X_GEO fails, try using Read10X
              x.mat <- Seurat::Read10X(data.dir = x)
              Seurat::CreateSeuratObject(counts = x.mat, project = basename(x))
            }, error = function(e) {
              tryCatch({
                # If Read10X fails, try using Read10X_h5
                x.mat <- Seurat::Read10X_h5(filename = x)
                if (is.list(x.mat)) {
                  # If result is a list, find the largest matrix
                  x.mat <- findLargestMatrix(x.mat)
                }
                Seurat::CreateSeuratObject(counts = x.mat, project = basename(x))
              }, error = function(e) {
                stop("Read10X_h5_Multi_Directory failed.")
              })
            })
          })
          seu.obj
        })
        
        if (isTRUE(merge)) {
          seu.obj <- mergeExperiments(seu.list)
        } else {
          seu.obj <- seu.list
        }
      } else if (!is.null(pf.count) && supp.type == "count") {
        seu.obj <- Seurat::CreateSeuratObject(counts = pf.count, project = acce)
      }
      return(seu.obj)
    }
  }
}

#' Extract Sample Metadata from GEO.
#'
#' @param acce GEO accession number.
#' @param platform Platform information/field. Default: NULL (all platforms).
#' @param ... Parameters for \code{\link{getGEO}}.
#'
#' @return Dataframe contains all metadata of provided GEO accession number.
#' @importFrom magrittr %>%
#' @importFrom GEOquery getGEO
#' @importFrom Biobase annotation experimentData pData phenoData notes sampleNames exprs
#' @export
#'
#' @examples
#' \donttest{
#' # users may need to set the size of the connection buffer
#' # Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 60)
#' # extract metadata of specified platform
#' GSE200257.meta <- ExtractGEOMeta(acce = "GSE200257", platform = "GPL24676")
#' }
ExtractGEOMeta <- function(acce, platform = NULL, ...) {
  # get GEO object
  if (is.null(platform)) {
    geo.obj <- GEOquery::getGEO(GEO = acce, ...)
    pfs <- sapply(geo.obj, function(x) {
      Biobase::annotation(x)
    })
    names(geo.obj) <- pfs
  } else {
    geo.obj <- list()
    pf.obj <- GEOobj(acce = acce, platform = platform, ...)
    geo.obj[[platform]] <- pf.obj
  }
  # extract metadata
  geo.meta.list <- lapply(names(geo.obj), function(x) {
    # extract general information
    pf.info <- ExtractGEOInfo(pf.obj = geo.obj[[x]], sample.wise = FALSE)
    pf.info$Platform <- x
    # select meta data
    pf.meta <- ExtractGEOSubMeta(pf.obj = geo.obj[[x]])
    pf.meta$Platform <- x
    # merge all dataframe
    pf.all <- merge(pf.meta, pf.info, by = "Platform", all.x = TRUE) %>% as.data.frame()
    pf.all
  })
  # all metadta
  if (length(geo.meta.list) > 1) {
    geo.meta.df <- data.table::rbindlist(geo.meta.list, fill = TRUE) %>% as.data.frame()
  } else {
    geo.meta.df <- geo.meta.list[[1]]
  }
  return(geo.meta.df)
}

# connect to GEO, extract GEO object, extract platform object
GEOobj <- function(acce, platform, ...) {
  # obtain GEO object
  geo.obj <- GEOquery::getGEO(GEO = acce, ...)

  # extract platform
  pfs <- sapply(geo.obj, function(x) {
    Biobase::annotation(x)
  })
  if (!platform %in% pfs) {
    stop(paste("The platform you provides is not valid!", paste(pfs, collapse = ", ")))
  }
  # extract platform data
  pf.idx <- which(pfs == platform[1])
  pf.obj <- geo.obj[[pf.idx]]
  return(pf.obj)
}

# merge cols into one
SimpleCol <- function(df, col) {
  cols <- grep(pattern = col, x = colnames(df), value = T)
  df.cols <- df[cols]
  value <- paste(c(t(df.cols)), collapse = ". ")
  return(value)
}

#' Extract GEO Study Information.
#'
#' @param pf.obj GEO object of platform.
#' @param sample.wise Logical value, whether to extract sample-wise information. Default: FALSE.
#'
#' @return A dataframe.
#'
ExtractGEOInfo <- function(pf.obj, sample.wise = FALSE) {
  # platform information
  pf.info <- Biobase::experimentData(pf.obj)
  # additional information
  pf.info.add.df <- as.data.frame(Biobase::pData(Biobase::phenoData(pf.obj)))
  if (sample.wise) {
    pf.info.final <- cbind(data.frame(
      Title = pf.info@title,
      Type = Biobase::notes(pf.info)$type,
      Abstract = pf.info@abstract,
      Design = pf.info@other$overall_design,
      SampleCount = length(Biobase::sampleNames(pf.obj)),
      SupplementaryFile = gsub(pattern = "\n", replacement = ", ", pf.info@other$supplementary_file),
      PMID = gsub(pattern = "\n", replacement = ", ", pf.info@pubMedIds), pf.info.add.df
    )) %>%
      as.data.frame()
  } else {
    used.cols <- c("organism", "molecule", "strategy", "extract_protocol", "data_processing")
    pf.info.add.used <- pf.info.add.df[, grep(pattern = paste(used.cols, collapse = "|"), colnames(pf.info.add.df), value = T)]
    pf.info.add.sim <- apply(pf.info.add.used, 2, function(x) {
      paste(unique(x), collapse = ", ")
    }) %>%
      t() %>%
      as.data.frame()
    # final information
    pf.info.final <- data.frame(
      Title = pf.info@title,
      Type = Biobase::notes(pf.info)$type,
      Organism = SimpleCol(df = pf.info.add.sim, col = "organism"),
      Abstract = pf.info@abstract,
      Design = pf.info@other$overall_design,
      SampleCount = length(Biobase::sampleNames(pf.obj)),
      Molecule = SimpleCol(df = pf.info.add.sim, col = "molecule"),
      ExtractProtocol = SimpleCol(df = pf.info.add.sim, col = "extract_protocol"),
      LibraryStrategy = SimpleCol(df = pf.info.add.sim, col = "strategy"),
      DataProcessing = SimpleCol(df = pf.info.add.sim, col = "data_processing"),
      SupplementaryFile = gsub(pattern = "\n", replacement = ", ", pf.info@other$supplementary_file),
      Contact = paste(pf.info@name, pf.info@contact, sep = "; "),
      PMID = gsub(pattern = "\n", replacement = ", ", pf.info@pubMedIds)
    )
  }
  return(pf.info.final)
}

#' Extract Sample Metadata.
#'
#' @param pf.obj GEO object of platform.
#'
#' @return A dataframe.
#'
ExtractGEOSubMeta <- function(pf.obj) {
  # extract sample detail information
  pf.info <- as.data.frame(Biobase::pData(Biobase::phenoData(pf.obj)))
  # select used basic cols
  valid.cols <- intersect(colnames(pf.info), c(c("title", "geo_accession", "source_name_ch1", "description")))
  pf.info.used <- pf.info[valid.cols]
  # process characteristics
  pf.info.charac <- pf.info[grep(pattern = "^characteristics", x = colnames(pf.info))]
  ## modify colnames
  pf.info.charac.colnames <-
    unique(apply(pf.info.charac, 2, function(x) {
      gsub(pattern = "(.*?): (.*)", replacement = "\\1", x = x)
    }))
  colnames(pf.info.charac) <- pf.info.charac.colnames
  ## modify values
  pf.info.charac <- apply(pf.info.charac, 2, function(x) {
    gsub(pattern = "(.*?): (.*)", replacement = "\\2", x = x)
  })
  ## final meta
  pf.meta <- cbind(pf.info.used, pf.info.charac) %>% as.data.frame()

  return(pf.meta)
}


#' Extract Raw Count Matrix from Supplementary Files.
#'
#' @param acce GEO accession number.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#'
#' @return A dataframe.
#'
ExtractGEOExpSupp <- function(acce, timeout = 3600, supp.idx = 1,manual.download= NULL) {
  # create tmp folder
  tmp.folder <- tempdir()
  # get current timeout
  if (!is.null(timeout)) {
    message("Change Timeout to: ", timeout)
    env.timeout <- getOption("timeout")
    on.exit(options(timeout = env.timeout)) # restore timeout
    options(timeout = timeout)
  }
  
  clean_directory_if_not_empty(file.path(tmp.folder, acce))
  
  # download supplementary file
  # supp.down.log <- GEOquery::getGEOSuppFiles(GEO = acce, baseDir = tmp.folder)
  if (!is.null(manual.download)) {
    if (!dir.exists(file.path(tmp.folder, acce))) {
      dir.create(path = file.path(tmp.folder, acce), recursive = TRUE)
    }
    file.copy(file.path(manual.download,"/"), file.path(tmp.folder, acce),overwrite = TRUE,recursive = T)
    raw_tar_files <- list.files(file.path(tmp.folder, acce), full.names = TRUE, pattern = ".tar",recursive = T)
  }else{
    supp.down.log <- tryCatch(
      expr = {
        GEOquery::getGEOSuppFiles(GEO = acce, baseDir = tmp.folder)
      },
      error = function(e) {
        print(e)
        stop("You can change the timeout with a larger value.")
      }
    )
    raw_tar_files <- rownames(supp.down.log)[grepl(".tar", rownames(supp.down.log))]
  }
  
  
  
  
  for (file_path in raw_tar_files) {
    # if(grepl(".tar.gz$",file_path)){
    #   gunzip(file_path,overwrite  = T)
    #   file_path <- gsub(pattern = "\\.gz", replacement = "", x = file_path)
    # }
    if(grepl(".tar",file_path)){
      extractTar(file_path)
    }
  }
  
  
  supp.file.paths = list.files(file.path(tmp.folder, acce), full.names = TRUE,recursive = T)
  supp.idx = supp.file.paths %>% sapply(isCountsMatrix)
  
  log_info(str_glue("{supp.file.paths}: {supp.idx}"))
  
  supp.file.paths = list.files(file.path(tmp.folder, acce), full.names = TRUE,recursive = T)
  if (sum(supp.idx) == 0) {
    stop("The counts files do not meet the criteria.")
  }else if (sum(supp.idx) == 1) {
    supp.file.path <- supp.file.paths[supp.idx]
  }else{
    supp.file.paths <- supp.file.paths[supp.idx]
    supp.file.path = file.path(tmp.folder,acce,"count.tar")
    tryCatch({
      tar(tarfile = supp.file.path, files = supp.file.paths)
      log_info("Files {str_c(supp.file.paths, collapse = ', ')} have been compressed into {supp.file.path}.")
    }, warning = function(w) {
      log_warn("Warning during compression: {w$message}")
      warning(w)
    }, error = function(e) {
      log_error("Error during compression: {e$message}")
      stop(e)
    })
  }
  supp.file.name = basename(supp.file.path)
  log_info(str_glue("reading counts from {supp.file.name}"))
  # remove unused supplementary file
  unused.supp <- setdiff(supp.file.paths, supp.file.path)
  unused.supp.remove <- file.remove(unused.supp)
  # file unzip
  file.ext <- tools::file_ext(supp.file.path)
  if (file.ext == "gz") {
    # gunzip file
    GEOquery::gunzip(supp.file.path, overwrite = TRUE)
    supp.file.path <- gsub(pattern = "\\.gz", replacement = "", x = supp.file.path)
    # update file extension
    file.ext <- tools::file_ext(supp.file.path)
  }
  if (file.ext == "tar"){
    # untar
    extractTar(supp.file.path,exdir =file.path(tmp.folder,acce,"sample") )
    # unzip
    unzip.log <- sapply(
      list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE,recursive = TRUE, pattern = "gz$"),
      function(x) {
        GEOquery::gunzip(x, overwrite = TRUE)
      }
    )
    # read files
    count.list <- lapply(
      list.files(file.path(tmp.folder, acce, "sample"),full.names = TRUE,recursive = TRUE),
      function(x) {
        sample.count <- data.table::fread(file = x) %>% as.data.frame()
        colnames(sample.count) <- c("GeneName", 
                                    paste0(
                                      gsub(pattern = "(GSM[0-9]*).*", replacement = "\\1", x = basename(x)),
                                      "__",
                                      colnames(sample.count)[-1]
                                    )
                                  )
        sample.count
      }
    )
    # create count matrix
    count.mat <- Reduce(f = function(x, y) {
      merge.mat <- merge(x, y, by = "GeneName", all = T)
    }, x = count.list)
    rownames(count.mat) <- count.mat$GeneName
    count.mat$GeneName <- NULL
  }else {
    if (file.ext %in% c("xlsx", "xls")) {
      # read excel file
      count.mat <- openxlsx::read.xlsx(xlsxFile = supp.file.path, rowNames = TRUE)
    } else if (file.ext %in% c("csv", "tsv", "txt")) {
      # read text file
      count.mat <- data.table::fread(file = supp.file.path) %>% as.data.frame()
      # the first column must be gene
      rownames(count.mat) <- count.mat[, 1]
      count.mat[, 1] <- NULL
    }
  }
  
  
  
  return(count.mat)
}


#' Fortmat Supplementary Files to 10x.
#'
#' @param acce GEO accession number.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param out.folder Output folder to save 10x files. Default: NULL (current working directory).
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}.
#' Default: TURE.
#'
#' @return NULL
#'
ExtractGEOExpSupp10x <- function(acce, supp.idx = 1, timeout = 3600,manual.download=NULL,
                                 out.folder = NULL, gene2feature = TRUE) {
  # create tmp folder
  tmp.folder <- tempdir()
  # get current timeout
  if (!is.null(timeout)) {
    message("Change Timeout to: ", timeout)
    env.timeout <- getOption("timeout")
    on.exit(options(timeout = env.timeout)) # restore timeout
    options(timeout = timeout)
  }
  # download supp file
  if (!is.null(manual.download)) {
    tryCatch({
      tar(tarfile = file.path(tmp.folder, acce,"manual.tar.gz"), files = manual.download)
      log_info("Files {str_c(supp.file.paths, collapse = ', ')} have been compressed into {supp.file.path}.")
    }, warning = function(w) {
      log_warn("Warning during compression: {w$message}")
      warning(w)
    }, error = function(e) {
      log_error("Error during compression: {e$message}")
      stop(e)
    })
    all.files = list.files(file.path(tmp.folder, acce), full.names = TRUE,recursive = T)
  }else{
    supp.down.log <- tryCatch(
      expr = {
        GEOquery::getGEOSuppFiles(GEO = acce, baseDir = tmp.folder)
      },
      error = function(e) {
        print(e)
        stop("You can change the timeout with a larger value.")
      }
    )
    all.files = row.names(supp.down.log)
    # check supp.idx
    if (supp.idx > nrow(supp.down.log)) {
      stop("Please provide valid supplementary file index.")
    }
  }
  # get used supplementary file
  supp.file.path <- all.files[supp.idx]
  # remove unused supplementary file
  unused.supp <- setdiff(all.files, supp.file.path)
  unused.supp.remove <- file.remove(unused.supp)
  # file unzip
  file.ext <- tools::file_ext(supp.file.path)
  if (file.ext == "gz") {
    # gunzip file
    GEOquery::gunzip(supp.file.path, overwrite = TRUE)
    supp.file.path <- gsub(pattern = "\\.gz", replacement = "", x = supp.file.path)
    # update file extension
    file.ext <- tools::file_ext(supp.file.path)
  }
  if (file.ext == "tar") {
    # untar
    utils::untar(supp.file.path, exdir = file.path(tmp.folder, acce, "sample"))
    # recognize valid files: barcodes.tsv.gz, genes.tsv.gz, matrix.mtx.gz and features.tsv.gz
    valid.pat <- "barcodes.tsv.gz$|genes.tsv.gz$|matrix.mtx.gz$|features.tsv.gz$|h5$"
    all.files <- list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = valid.pat)
    # change file name
    if (gene2feature) {
      change.name.log <- sapply(all.files, function(x) {
        if (grepl(pattern = "genes.tsv.gz$", x = x)) {
          new.name <- gsub(pattern = "genes.tsv.gz$", replacement = "features.tsv.gz", x = x)
          file.rename(from = x, to = new.name)
        }
      })
      all.files <- list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = valid.pat)
    }
    # prepare out folder
    if (is.null(out.folder)) {
      out.folder <- getwd()
    }
    # get folder
    all.sample.folder <- sapply(all.files, function(x) {
      # get basename and dirname
      file.name <- basename(x)
      dir.name <- dirname(x)
      # remove file type tag
      file.name <- gsub(pattern = valid.pat, replacement = "", x = file.name)
      # remove possible _ and .
      file.name <- gsub(pattern = "[_.]$", replacement = "", x = file.name)
      file.folder <- file.path(out.folder, file.name)
    })
    # create folder and move file
    move.file.log <- sapply(all.files, function(x) {
      # get folder name
      folder.name <- all.sample.folder[x]
      # create folder
      if (!dir.exists(folder.name)) {
        dir.create(path = folder.name, recursive = TRUE)
      }
      new.file.name <- gsub(pattern = paste0(".*(barcodes.tsv.gz$|genes.tsv.gz$|matrix.mtx.gz$|features.tsv.gz$|h5$)"), replacement = "\\1", x = x)
      # move file
      copy.tag <- file.copy(from = x, to = file.path(folder.name, new.file.name))
      # remove the original file
      remove.tag <- file.remove(x)
      copy.tag
    })
    message("Process 10x fiels done! All files are in ", out.folder)
  } else {
    log_warn("Detect non-tar file for 10x mode!")
    log_info("try 10x single mode")
    
    tryCatch(
      expr = {
        ExtractGEOExpSupp10xSingle(acce = acce, timeout = timeout,manual.download= manual.download, out.folder = out.folder, gene2feature = gene2feature)
      },
      error = function(e) {
        print(e)
        stop("fail to try 10x single mode!")
      }
    )
  }
}

#' Fortmat Supplementary Files to 10x (separate files).
#'
#' @param acce GEO accession number.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param out.folder Output folder to save 10x files. Default: NULL (current working directory).
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}.
#' Default: TURE.
#'
#' @return NULL
#'
ExtractGEOExpSupp10xSingle <- function(acce, timeout = 3600, out.folder = NULL, gene2feature = TRUE,manual.download=NULL) {
  # create tmp folder
  tmp.folder <- tempdir()
  # get current timeout
  if (!is.null(timeout)) {
    message("Change Timeout to: ", timeout)
    env.timeout <- getOption("timeout")
    on.exit(options(timeout = env.timeout)) # restore timeout
    options(timeout = timeout)
  }
  
  if (!(is.null(manual.download))) {
    file.copy(file.path(manual.download,"/"), file.path(tmp.folder, acce,"/"),overwrite = TRUE,recursive = T)
  }else{
  # download supp file
    supp.down.log <- tryCatch(
      expr = {
        GEOquery::getGEOSuppFiles(GEO = acce, baseDir = tmp.folder)
      },
      error = function(e) {
        print(e)
        stop("You can change the timeout with a larger value.")
      }
    )
  }
  # get valid files
  valid.pat <- "barcodes.tsv.gz$|genes.tsv.gz$|mtx.gz$|features.tsv.gz$"
  all.files <- list.files(file.path(tmp.folder, acce), full.names = TRUE, pattern = valid.pat,recursive = T)
  # change file name
  if (gene2feature) {
    change.name.log <- sapply(all.files, function(x) {
      if (grepl(pattern = "genes.tsv.gz$", x = x)) {
        new.name <- gsub(pattern = "genes.tsv.gz$", replacement = "features.tsv.gz", x = x)
        file.rename(from = x, to = new.name)
      }
    })
    all.files <- list.files(file.path(tmp.folder, acce), full.names = TRUE, pattern = valid.pat)
  }
  # prepare out folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  
  if (length(all.files)!=3) {
    log_warn(str_glue("10X single should be 3 files"))
    all.files_base <- list.files(file.path(tmp.folder, acce), full.names = F, pattern = valid.pat)
    log_warn(str_glue("Files including: {str_c(all.files_base, collapse = '\n ')}."))
    if(length(all.files)>3) {
      
      log_warn(str_glue("Automatically Choosing files..."))
      ## only choosing mtx file
      valid.pat_mtx <-"mtx.gz"
      all.files_mtx = list.files(file.path(tmp.folder, acce), full.names = TRUE, pattern = valid.pat_mtx)
      
      for (file in all.files_mtx) {
        if (!isCountMTX(file)) {
          file.remove(file)
        }
      }
      all.files_mtx = list.files(file.path(tmp.folder, acce), full.names = TRUE, pattern = valid.pat_mtx)
      all.files <- list.files(file.path(tmp.folder, acce), full.names = TRUE, pattern = valid.pat)
    }
    
    if(length(all.files)>3) {
      all.files_mtx = list.files(file.path(tmp.folder, acce), full.names = TRUE, pattern = valid.pat_mtx)[1]
      log_warnings(str_glue("Choosing frist mtx:{all.files_mtx[1]}..."))
    }
      
    
    if (length(all.files)<3) {
      stop(str_glue("Files only including: {str_c(all.files_base, collapse = '\n ')}."))
    }
    file.rename(from = all.files_mtx,
                to = gsub(valid.pat_mtx,replacement = "matrix.mtx.gz",x = all.files_mtx))
    
    all.files <- list.files(file.path(tmp.folder, acce), full.names = TRUE, pattern = valid.pat)
    
  }
  
  # get folder
  all.sample.folder <- sapply(all.files, function(x) {
    # get basename and dirname
    file.name <- basename(x)
    dir.name <- dirname(x)
    # remove file type tag
    file.name <- gsub(pattern = valid.pat, replacement = "", x = file.name)
    # remove possible _ and .
    file.name <- gsub(pattern = "[_.]$", replacement = "", x = file.name)
    file.folder <- file.path(out.folder, file.name)
  })
  
  all.sample.folder = findCommonPart(all.sample.folder)
  # create folder and move file
  move.file.log <- sapply(all.files, function(x) {
    # get folder name
    folder.name <- all.sample.folder[x]
    # create folder
    if (!dir.exists(folder.name)) {
      dir.create(path = folder.name, recursive = TRUE)
    }
    new.file.name <- gsub(pattern = paste0(".*(barcodes.tsv.gz$|genes.tsv.gz$|matrix.mtx.gz$|features.tsv.gz$)"), replacement = "\\1", x = x)
    # move file
    copy.tag <- file.copy(from = x, to = file.path(folder.name, new.file.name),overwrite = T)
    # remove the original file
    # remove.tag <- file.remove(x)
    copy.tag
  })
  message("Process 10x fiels done! All files are in ", out.folder)
}


#' Extract Raw Count Matrix from Supplementary Files.
#'
#' @param acce GEO accession number.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#'
#' @return A dataframe.
#'
ExtractGEOExpSuppObj <- function(acce, timeout = 3600000, supp.idx = 1,manual.download = NULL,
                                 supp.type=c("Seurat","SingleCellExperiment","Loom","AnnData")) {
  # create tmp folder
  tmp.folder <- tempdir()
  # get current timeout
  if (!is.null(timeout)) {
    message("Change Timeout to: ", timeout)
    env.timeout <- getOption("timeout")
    on.exit(options(timeout = env.timeout)) # restore timeout
    options(timeout = timeout)
  }
  
  if (!is.null(manual.download)) {
    if (!dir.exists(file.path(tmp.folder, acce))) {
      dir.create(path = file.path(tmp.folder, acce), recursive = TRUE)
    }
    file.copy(file.path(manual.download,"/"), file.path(tmp.folder, acce),overwrite = TRUE,recursive = T)
    all.files <- list.files(file.path(tmp.folder, acce), full.names = TRUE,recursive = T)
   
  }else{
    # download supplementary file
    # supp.down.log <- GEOquery::getGEOSuppFiles(GEO = acce, baseDir = tmp.folder)
    supp.down.log <- tryCatch(
      expr = {
        GEOquery::getGEOSuppFiles(GEO = acce, baseDir = tmp.folder)
      },
      error = function(e) {
        print(e)
        stop("You can change the timeout with a larger value.")
      }
    )
    if (supp.idx > nrow(supp.down.log)) {
      stop("Please provide valid supplementary file index.")
    }
    all.files = rownames(supp.down.log)
  }
  supp.file.path <- all.files[supp.idx]
  isValidFile <- checkFileExtension(supp.file.path)
  
  # Output the result using logger
  if (isValidFile) {
    log_info("Has provided a valid supplementary file index")
  } else {
    log_warnings("Please provide valid supplementary file index.Try recongize auto ...")
    
    supp.file.path =  findOBJFile(supp.type, all.files)
    
    isValidFile <- checkFileExtension(supp.file.path)
    
    if (!isValidFile) {
      log_error(str_glue("Cannot find {supp.type} file,please provide valid supplementary file index or check files. "))
    }
  }
  
  # remove unused supplementary file
  unused.supp <- setdiff(all.files, supp.file.path)
  # unused.supp.remove <- file.remove(unused.supp)
  if(length(supp.file.path)> 1 ){
    stop(str_glue("More than one {supp.type} files, waiting update pipeline."))
  }
  # file unzip
  file.ext <- tools::file_ext(supp.file.path)
  if (file.ext == "gz") {
    # gunzip file
    GEOquery::gunzip(supp.file.path, overwrite = TRUE)
    supp.file.path <- gsub(pattern = "\\.gz", replacement = "", x = supp.file.path)
    # update file extension
    file.ext <- tools::file_ext(supp.file.path)
  }
  
  if (file.ext == "tar") {
    # untar
    utils::untar(supp.file.path, exdir = file.path(tmp.folder, acce, "sample"))
    # unzip
    unzip.log <- sapply(
      list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = "gz$"),
      function(x) {
        GEOquery::gunzip(x, overwrite = TRUE)
      }
    )
    
    supp.file.path <- tools::file_path_sans_ext(supp.file.path)
    
  } 
  
    if (tolower(file.ext) %in%  c("rds", "rdata", "h5ad", "loom")) {
      # read excel file
      # supp.type=c("Seurat","SingleCellExperiment","Loom","AnnData")
      
        outs_path = file.path(".",acce,"outs/processed.h5ad")
        
        if(supp.type == "Seurat"){
          seurat = readRDS(supp.file.path)
          
          seurat %>% ExportSeurat(to = "AnnData",conda.path = "/opt/conda",anndata.file = outs_path)
          
        }else if(supp.type == "SingleCellExperiment"){
          sce = readRDS(supp.file.path)
          
          sce %>% SCEAnnData(from = "SingleCellExperiment",to = "AnnData",conda.path = "/opt/conda",anndata.file = outs_path)
          
        }else if(supp.type == "Loom"){
          stop("Convert loom to Anndata fail!")
        }else if(supp.type == "AnnData"){
          file.rename(supp.file.path,outs_path)
        }
        
        log_trace(str_glue("Saved anndata in {outs_path}"))

      } else {

      stop(str_glue("Cannot find {supp.type} file,please provide valid supplementary file index or check files."))
      
    }
  return()
}



#' Extract Raw Count Matrix from Supplementary Files or Fortmat Supplementary Files to 10x.
#'
#' @param acce GEO accession number.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param supp.type The type of downloaded supplementary files, choose from count (count matrix file or single count matrix file),
#' 10x (cellranger output files in tar/gz supplementary files, contains barcodes, genes/features and matrix, e.g. GSE200257)
#' and 10xSingle (cellranger output files in supplementary files directly, e.g. GSE236082). Default: count.
#' @param out.folder Output folder to save 10x files. Default: NULL (current working directory).
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}. Default: TRUE.
#' Default: TURE.
#'
#' @return Count matrix (\code{supp.type} is count) or NULL (\code{supp.type} is 10x).
#'
ExtractGEOExpSuppAll <- function(acce, supp.idx = 1, timeout = 3600,manual.download= manual.download,
                                 supp.type = c("count", "10x", "10xSingle"), out.folder = NULL, gene2feature = TRUE) {
  if (supp.type == "count") {
    count.mat <- ExtractGEOExpSupp(acce = acce, supp.idx = supp.idx, manual.download= manual.download, timeout = timeout)
    return(count.mat)
  } else if (supp.type == "10x") {
    ExtractGEOExpSupp10x(acce = acce, supp.idx = supp.idx, timeout = timeout,manual.download= manual.download, out.folder = out.folder, gene2feature = gene2feature)
    return(NULL)
  } else if (supp.type == "10xSingle") {
    ExtractGEOExpSupp10xSingle(acce = acce, timeout = timeout,manual.download= manual.download, out.folder = out.folder, gene2feature = gene2feature)
    return(NULL)
  }
}

#' Extract Raw Count Matrix or Fortmat Supplementary Files to 10x.
#'
#' @param pf.obj GEO object of platform.
#' @param acce GEO accession number.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#' @param down.supp Logical value, whether to download supplementary files to create count matrix. If TRUE, always
#' download supplementary files. If FALSE, use \code{ExpressionSet} (If contains non-integer or emoty,
#' download supplementary files automatically). Default: FALSE.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param supp.type The type of downloaded supplementary files, choose from count (count matrix file or single count matrix file),
#' 10x (cellranger output files in tar/gz supplementary files, contains barcodes, genes/features and matrix, e.g. GSE200257)
#' and 10xSingle (cellranger output files in supplementary files directly, e.g. GSE236082). Default: count.
#' @param out.folder Output folder to save 10x files. Default: NULL (current working directory).
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}. Default: TRUE.
#'
#' @return Count matrix (\code{supp.type} is count) or NULL (\code{supp.type} is 10x/10xSingle).
#'
ExtractGEOExp <- function(pf.obj, acce, supp.idx = 1, down.supp = FALSE, timeout = 3600,manual.download= manual.download,
                          supp.type = c("count", "10x", "10xSingle"), out.folder = NULL, gene2feature = TRUE) {
  # check parameters
  supp.type <- match.arg(arg = supp.type)
  # download supplementary files
  if (down.supp) {
    exp.data <- ExtractGEOExpSuppAll(
      acce = acce, supp.idx = supp.idx, timeout = timeout,manual.download= manual.download,
      supp.type = supp.type, out.folder = out.folder, gene2feature = gene2feature
    )
  } else {
    expr.mat <- Biobase::exprs(pf.obj)
    if (nrow(expr.mat) == 0) {
      message("Matrix not available! Downloading supplementary files.")
      exp.data <- ExtractGEOExpSuppAll(
        acce = acce, supp.idx = supp.idx, timeout = timeout,manual.download= manual.download,
        supp.type = supp.type, out.folder = out.folder, gene2feature = gene2feature
      )
    } else {
      if (all(expr.mat %% 1 == 0)) {
        exp.data <- expr.mat
      } else {
        message("Matrix contains non-integer values! Downloading supplementary files.")
        exp.data <- ExtractGEOExpSuppAll(
          acce = acce, supp.idx = supp.idx, timeout = timeout,manual.download= manual.download,
          supp.type = supp.type, out.folder = out.folder, gene2feature = gene2feature
        )
      }
    }
  }
  return(exp.data)
}

ExtractGEOMetaFile = function(acce,workspace,meta_file_pattern){
  dir_meta = file.path(workspace ,"meta")
  # dir_meta = file.path(".") ## /GSE161382/GSE161382_metadata.txt.gz
  if (!dir.exists(dir_meta)) {
    dir.create(path = dir_meta, recursive = TRUE)
  }
  supp.down.log = GEOquery::getGEOSuppFiles(GEO = acce, baseDir = dir_meta,filter_regex = meta_file_pattern)
  list_meta = rownames(supp.down.log) %>% map(function(supp.file.path){
    # supp.file.path =   rownames(supp.down.log)[1]
    file.ext <- tools::file_ext(supp.file.path)
    if (file.ext == "gz") {
      # gunzip file
      GEOquery::gunzip(supp.file.path, overwrite = TRUE)
      supp.file.path <- gsub(pattern = "\\.gz", replacement = "", x = supp.file.path)
      # update file extension
      file.ext <- tools::file_ext(supp.file.path)
    }
    if (file.ext == "tar") {
      # untar
      utils::untar(supp.file.path, exdir = file.path(tmp.folder, acce, "sample"))
      # unzip
      unzip.log <- sapply(
        list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = "gz$"),
        function(x) {
          GEOquery::gunzip(x, overwrite = TRUE)
        }
      )
      # read files
      count.list <- lapply(
        list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE),
        function(x) {
          sample.count <- data.table::fread(file = x) %>% as.data.frame()
          colnames(sample.count) <- c("GeneName", gsub(pattern = "(GSM[0-9]*).*", replacement = "\\1", x = basename(x)))
          sample.count
        }
      )
      # create count matrix
      df_meta <- Reduce(f = function(x, y) {
        merge.mat <- merge(x, y, by = "GeneName", all = T)
      }, x = count.list)
      rownames(df_meta) <- df_meta$GeneName
      df_meta$GeneName <- NULL
    } else {
      if (file.ext %in% c("xlsx", "xls")) {
        # read excel file
        df_meta <- openxlsx::read.xlsx(xlsxFile = supp.file.path, rowNames = TRUE)
      } else if (file.ext %in% c("csv", "tsv", "txt")) {
        # read text file
        df_meta <- data.table::fread(file = supp.file.path) %>% as.data.frame()
        head(df_meta) %>% log_info()
        # the first column must be gene
        rownames(df_meta) <- df_meta[, 1]
        df_meta[, 1] <- NULL
      }
    }
  return(df_meta)
  })
  meta = list_meta %>% reduce(rbind)
  # log_info("Print head of meta data table...")
  # log_info(meta %>% head)
  return(meta)
}

