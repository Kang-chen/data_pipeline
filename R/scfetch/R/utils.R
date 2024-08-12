# merge Seurat object, modify from: https://github.com/dosorio/rPanglaoDB/blob/master/R/mergeExperiments.R
# library(jsonlite)
library(data.table)
library(readxl)
library(R.utils) # For gunzip
library(utils) # For untar

mergeExperiments <- function(experimentList) {
  el.df <- lapply(experimentList, FUN = function(x) {
    dim(x)
  })
  el.df <- as.data.frame(t(as.data.frame(el.df)))
  rownames(el.df) <- 1:nrow(el.df)
  el.df.vec <- apply(el.df, 1, function(row) all(row != 0))
  empty.seu.index <- names(el.df.vec)[el.df.vec == FALSE]
  if (length(empty.seu.index) > 0) {
    message("Detect empty SeuratObject: ", paste0(empty.seu.index, collapse = ", "), ". Skip these!")
  }
  experimentList <- experimentList[el.df.vec]
  for (i in seq_along(experimentList)[-1]) {
    experimentList[[1]] <- suppressWarnings(merge(experimentList[[1]], experimentList[[i]]))
    experimentList[[i]] <- methods::new("Seurat")
  }
  experimentList <- experimentList[[1]]
  return(experimentList)
}

# check columns existence
CheckColumns <- function(df, columns) {
  if (!all(columns %in% colnames(df))) {
    miss.cols <- setdiff(columns, colnames(df))
    stop(paste0(paste(miss.cols, collapse = ", "), " does not exist, Please Check!"))
  }
}

# used in UCSCCellBrowser and cellxgene, merge multiple attributes
PasteAttr <- function(df, attr) {
  for (at in attr) {
    df[[at]] <- sapply(df[[at]], function(x) {
      paste0(x, collapse = ", ")
    })
  }
  return(df)
}

# used in UCSCCellBrowser, recursively extract samples
ExtractSample <- function(df, base.url, json.folder, quiet) {
  # prepare json
  if (base.url != json.folder) {
    df.json <- file.path(base.url, df$name, "dataset.json")
    names(df.json) <- df$name
    df.json.folder <- file.path(json.folder, df$name)
    names(df.json.folder) <- df$name
    df.desc <- file.path(base.url, df$name, "desc.json")
    names(df.desc) <- df$name
    dird <- sapply(df.json.folder, function(x) {
      dir.create(x, showWarnings = FALSE, recursive = TRUE)
    })
    down.status <- lapply(df$name, function(x) {
      utils::download.file(url = df.json[x], destfile = file.path(df.json.folder[x], "dataset.json"), quiet = quiet, mode = "wb", method = "wget", extra = "--no-check-certificate")
      utils::download.file(url = df.desc[x], destfile = file.path(df.json.folder[x], "desc.json"), quiet = quiet, mode = "wb", method = "wget", extra = "--no-check-certificate")
    })
  }
  # process
  if (!"isCollection" %in% colnames(df)) {
    return(df)
  } else {
    cf <- df[!is.na(df$isCollection), ]
    sf <- df[is.na(df$isCollection), ]
    cu.json.folder <- file.path(json.folder, cf$name)
    cul <- lapply(file.path(cu.json.folder, "dataset.json"), function(x) {
      x.json <- jsonlite::fromJSON(txt = x)
      x.df <- jsonlite::flatten(x.json$datasets)
      colnames(x.df) <- gsub(pattern = ".*\\.", replacement = "", x = colnames(x.df))
      x.df
    })
    cu.df <- data.table::rbindlist(cul, fill = TRUE)
    # df = data.table::rbindlist(list(sf, cu.df), fill = TRUE)
    # return(list(sf, ExtractSample(cu.df)))
    return(data.table::rbindlist(list(sf, ExtractSample(df = cu.df, base.url = base.url, json.folder = json.folder, quiet = quiet)), fill = TRUE))
  }
}

# used in UCSCCellBrowser, recursively extract samples online
ExtractSampleOnline <- function(df) {
  base.url <- "https://cells.ucsc.edu/"
  if (!"isCollection" %in% colnames(df)) {
    return(df)
  } else {
    cf <- df[!is.na(df$isCollection), ]
    sf <- df[is.na(df$isCollection), ]
    cu <- file.path(base.url, cf$name, "dataset.json")
    cul <- lapply(cu, function(x) {
      x.json <- jsonlite::fromJSON(txt = x)
      x.df <- jsonlite::flatten(x.json$datasets)
      colnames(x.df) <- gsub(pattern = ".*\\.", replacement = "", x = colnames(x.df))
      x.df
    })
    cu.df <- data.table::rbindlist(cul, fill = TRUE)
    # df = data.table::rbindlist(list(sf, cu.df), fill = TRUE)
    # return(list(sf, ExtractSample(cu.df)))
    return(data.table::rbindlist(list(sf, ExtractSampleOnline(cu.df)), fill = TRUE))
  }
}

# used in UCSCCellBrowser, inherit attributes from parents
InheritParient <- function(df, attr) {
  for (at in attr) {
    df[[at]] <- ifelse(df[[at]] == "", ifelse(df[[paste0("parent_", at)]] == "", "", paste0(df[[paste0("parent_", at)]], "|parent")), df[[at]])
  }
  return(df)
}

# used in UCSCCellBrowser, extract sample attribute
ExtractDesc <- function(lst, attr) {
  at.list <- list()
  for (atn in names(attr)) {
    at <- attr[atn]
    if (at %in% names(lst)) {
      at.value <- paste0(lst[[at]], collapse = ", ")
    } else {
      at.value <- ""
    }
    at.list[[atn]] <- at.value
  }
  return(as.data.frame(at.list))
}

# used in UCSCCellBrowser, filter attributes and return index
CheckParas <- function(df, column, para.value, fuzzy.match = TRUE) {
  # convert to lower case to avoid case-sensitive
  all.values <- gsub("\\|parent$", replacement = "", x = tolower(df[[column]]))
  if (is.null(para.value)) {
    message("Use all ", column, " as input!")
    # return row index
    value <- 1:nrow(df)
  } else {
    para.value <- tolower(para.value)
    # deal with fuzzy match
    if (fuzzy.match) {
      value.list <- sapply(para.value, function(x) grep(pattern = x, x = all.values, fixed = TRUE))
      value.list.len <- sapply(value.list, function(x) length(x))
      invalid.value <- names(value.list.len[value.list.len == 0])
      value <- unique(unlist(value.list))
    } else {
      if (column %in% c("body_parts", "diseases", "organisms", "projects")) {
        # value contains dot
        value.list <- list()
        for (pv in para.value) {
          pv.vec <- c()
          for (avi in 1:length(all.values)) {
            av <- all.values[avi]
            av.vec <- strsplit(x = av, split = ", ")[[1]]
            if (pv %in% av.vec) {
              pv.vec <- c(pv.vec, avi)
            }
          }
          value.list[[pv]] <- pv.vec
        }
        # get invalid value
        invalid.value <- setdiff(para.value, names(value.list))
        value <- unique(unlist(value.list))
      } else {
        # value doesn't contain dot
        value.list <- lapply(para.value, function(x) {
          which(all.values %in% x)
        })
        names(value.list) <- para.value
        value.list.len <- sapply(value.list, function(x) length(x))
        invalid.value <- names(value.list.len[value.list.len == 0])
        value <- unique(unlist(value.list))
      }
    }
    # print invalid value
    if (length(invalid.value) > 0) {
      message(paste0(invalid.value, collapse = ", "), " are not valid for ", column)
    }
    # deal with empty value
    if (length(value) == 0) {
      warning("There is no valid value under ", paste(para.value, collapse = ", "), " in ", column, ". This filter is invalid!")
      value <- 1:nrow(df)
    }
  }
  return(value)
}

# used in UCSCCellBrowser, create seurat object (add coord to metadata)
Load2Seurat <- function(exp.file, barcode.url = NULL, feature.url = NULL,
                        meta.file, coord.file = NULL, name = NULL, obs.value.filter = NULL,
                        obs.keys = NULL, include.genes = NULL) {
  # source: https://cellbrowser.readthedocs.io/en/master/load.html
  # read matrix
  if (is.null(barcode.url)) {
    mat <- data.table::fread(exp.file, check.names = FALSE)
    # get genes
    genes <- gsub(".+[|]", "", mat[, 1][[1]])
    # modify mat genes
    mat <- data.frame(mat[, -1], row.names = genes, check.names = FALSE)
  } else {
    # with default parameters
    mat <- Read10XOnline(matrix.url = exp.file, barcode.url = barcode.url, feature.url = feature.url)
  }
  # read metadata
  meta <- data.frame(data.table::fread(meta.file, check.names = FALSE), row.names = 1, check.names = FALSE)
  # filter dataset metadata
  ## filter cell's metadata values
  if (!is.null(obs.value.filter)) {
    meta <- meta %>% dplyr::filter(eval(rlang::parse_expr(obs.value.filter)))
  }
  ## filter cell's metadata colnames
  if (!is.null(obs.keys)) {
    meta <- meta[obs.keys]
  }
  ## filter genes
  if (!is.null(include.genes)) {
    mat <- mat[include.genes, rownames(meta)]
  } else {
    mat <- mat[, rownames(meta)]
  }
  if (is.null(coord.file)) {
    seu.obj <- Seurat::CreateSeuratObject(counts = mat, project = name, meta.data = meta)
  } else {
    # prepare coord file
    coord.list <- lapply(1:length(coord.file), function(x) {
      coord.name <- gsub(pattern = ".coords.tsv.gz", replacement = "", x = basename(coord.file[x]))
      coord.df <- data.frame(data.table::fread(coord.file[x], check.names = FALSE), row.names = 1, check.names = FALSE)
      colnames(coord.df) <- paste(coord.name, 1:ncol(coord.df), sep = "_")
      coord.df$Barcode <- rownames(coord.df)
      return(coord.df)
    })
    if (length(coord.file) == 1) {
      all.coord.df <- coord.list[[1]]
    } else {
      all.coord.df <- coord.list %>% purrr::reduce(dplyr::full_join, by = "Barcode")
    }
    # merge metadata
    meta <- merge(meta, all.coord.df, by.x = 0, by.y = "Barcode", all.x = TRUE) %>%
      tibble::column_to_rownames(var = "Row.names")
    seu.obj <- Seurat::CreateSeuratObject(counts = mat, project = name, meta.data = meta)
  }
  return(seu.obj)
}

Loading2DESeq2 <- function(mat, meta, fmu) {
  # loadding into DESeq2
  if (is.null(meta)) {
    meta <- data.frame(condition = colnames(mat))
    rownames(meta) <- colnames(mat)
    meta$condition <- as.factor(meta$condition)
    fmu <- "condition"
  } else {
    if (all(rownames(meta) != colnames(mat))) {
      stop("The columns of the count matrix and the rows of the meta.data are not in the same order!")
    }
  }
  if (is.null(fmu)) {
    message("The condition column (fmu) is empty, use the first column!")
    fmu <- colnames(meta)[1]
  }
  fmu.used <- stats::formula(paste("~", fmu))
  de.obj <- DESeq2::DESeqDataSetFromMatrix(
    countData = mat,
    colData = meta, design = fmu.used
  )
  return(de.obj)
}

# used in UCSCCellBrowser
# source: https://github.com/satijalab/seurat/blob/master/R/utilities.R#L1949
ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}

# used in UCSCCellBrowser, load cellranger output to matrix
Read10XOnline <- function(matrix.url, barcode.url, feature.url, gene.column = 2,
                          cell.column = 1, unique.features = TRUE, strip.suffix = FALSE) {
  # load matrix
  data <- Matrix::readMM(file = gzcon(url(matrix.url)))
  # load barcode
  cell.barcodes <- as.data.frame(data.table::fread(barcode.url, header = FALSE))
  cn <- ifelse(ncol(x = cell.barcodes) > 1, cell.column, 1)
  cell.names <- cell.barcodes[, cn]
  if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
    cell.names <- as.vector(x = as.character(x = sapply(
      X = cell.names,
      FUN = ExtractField,
      field = 1,
      delim = "-"
    )))
  }
  # matrix colnames
  colnames(x = data) <- cell.names
  # load feature
  feature.names <- as.data.frame(data.table::fread(feature.url, header = FALSE))
  # modify gene column
  gene.column <- min(ncol(feature.names), gene.column)
  if (any(is.na(x = feature.names[, gene.column]))) {
    warning(
      "Some features names are NA. Replacing NA names with ID from the opposite column requested",
      call. = FALSE,
      immediate. = TRUE
    )
    na.features <- which(x = is.na(x = feature.names[, gene.column]))
    replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
    feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
  }

  # modify matrix rownames
  if (unique.features) {
    rownames(x = data) <- make.unique(names = feature.names[, gene.column])
  } else {
    rownames(x = data) <- feature.names[, gene.column]
  }
  # In cell ranger 3.0, a third column specifying the type of data was added
  # and we will return each type of data as a separate matrix
  if (ncol(x = feature.names) > 2) {
    data_types <- factor(x = feature.names$V3)
    lvls <- levels(x = data_types)
    if (length(x = lvls) > 1) {
      message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
    }
    expr_name <- "Gene Expression"
    if (expr_name %in% lvls) { # Return Gene Expression first
      lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
    }
    data <- lapply(
      X = lvls,
      FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      }
    )
    names(x = data) <- lvls
  } else {
    data <- list(data)
  }
  # convert to dgCMatrix
  final.data <- Seurat::as.sparse(data[[1]])
  return(final.data)
}

# used in GEO, check the integrity of 10x files
Check10XFiles <- function(folders, gene2feature) {
  folders.flag <- sapply(folders, function(x) {
    if (gene2feature) {
      if (file.exists(file.path(x, "matrix.mtx.gz")) &&
        file.exists(file.path(x, "barcodes.tsv.gz")) &&
        file.exists(file.path(x, "features.tsv.gz"))) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      if (file.exists(file.path(x, "matrix.mtx.gz")) &&
        file.exists(file.path(x, "barcodes.tsv.gz"))) {
        if (file.exists(file.path(x, "features.tsv.gz")) ||
          file.exists(file.path(x, "genes.tsv.gz"))) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      } else if(file.exists(file.path(x, "h5"))){
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
  })
  valid.folders <- folders[folders.flag]
  drop.folders <- setdiff(folders, valid.folders)
  if (length(drop.folders) > 0) {
    if (gene2feature) {
      message(paste0(drop.folders, collapse = ", "), " don't contain matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz.")
    } else {
      message(paste0(drop.folders, collapse = ", "), " don't contain matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz/genes.tsv.gz.")
    }
  }
  return(valid.folders)
}

# check pair-end or single-end
CheckBam <- function(bam, samtools.path) {
  samtools.cmd <- paste(samtools.path, "view -h", bam, "2>/dev/null |head -n 100000 |", samtools.path, "view -c -f 1 -")
  # run command
  message(paste("Check pair/single end: ", samtools.cmd))
  samtools.status <- system(samtools.cmd, intern = TRUE)
  samtools.status <- as.numeric(samtools.status)
  if (samtools.status == 0) {
    return(FALSE)
  } else if (samtools.status > 0) {
    return(TRUE)
  }
}

# used in cellxxgene, extract content from url
URLRetrieval <- function(url) {
  url.page <- curl::curl_fetch_memory(url)
  url.content <- jsonlite::fromJSON(rawToChar(url.page$content))
  return(url.content)
}

# used in cellxgene, merge multiple attributes
PasteAttrCXG <- function(df, attr, col) {
  for (at in attr) {
    df[[at]] <- sapply(df[[at]], function(x) {
      paste0(x[[col]], collapse = ", ")
    })
  }
  return(df)
}

# used in cellxgene, filter datasets
cellxgeneAttrFilter <- function(df, attr, attr.value) {
  if (is.null(attr.value)) {
    message("Use all ", attr, " as input!")
    # return row index
    value <- 1:nrow(df)
  } else {
    # lower case
    attr.value <- tolower(attr.value)
    all.values <- df[[attr]]
    # value contains dot
    value.list <- list()
    for (pv in attr.value) {
      pv.vec <- c()
      for (avi in 1:length(all.values)) {
        av <- all.values[avi]
        av.vec <- strsplit(x = av, split = ", ")[[1]] %>% tolower()
        if (pv %in% av.vec) {
          pv.vec <- c(pv.vec, avi)
        }
      }
      value.list[[pv]] <- pv.vec
    }
    # get invalid value
    invalid.value <- setdiff(attr.value, names(value.list))
    value <- unique(unlist(value.list))
    # print invalid value
    if (length(invalid.value) > 0) {
      message(paste0(invalid.value, collapse = ", "), " are not valid for ", attr)
    }
    # deal with empty value
    if (length(value) == 0) {
      warning("There is no valid value under ", paste(attr.value, collapse = ", "), " in ", attr, ". This filter is invalid!")
      value <- 1:nrow(df)
    }
  }
  return(value)
}

# used in cellxgene, get download urls
PostDatasetURL <- function(url) {
  response <- httr::POST(url)
  httr::stop_for_status(response)
  result <- httr::content(response, as = "text", encoding = "UTF-8")
  result.list <- jsonlite::fromJSON(result)
  presigned.url <- result.list$presigned_url
  names(presigned.url) <- result.list$file_name
  return(presigned.url)
}

# used in cellxgene, prepare download urls with metadata
# reference: https://gist.github.com/ivirshup/f1a1603db69de3888eacb4bdb6a9317a
PrepareCELLxGENEUrls <- function(df, fe) {
  # file extension column
  fe.id <- paste0(fe, "_id")
  CheckColumns(df = df, columns = c("dataset_id", fe.id, "dataset_description"))
  invalid.df <- df[is.na(df[[fe.id]]) | is.na(df$dataset_id) | df$dataset_id == "" | df[[fe.id]] == "", ]
  message("Detect ", nrow(invalid.df), " invalid metadata (", fe.id, "/dataset_id is empty or NA).")
  valid.df <- df[!(is.na(df[[fe.id]]) | is.na(df$dataset_id) | df$dataset_id == "" | df[[fe.id]] == ""), ]
  valid.urls <- df[[fe.id]]
  valid.names <- make.names(valid.df$dataset_description, unique = TRUE)
  valid.filenames <- paste0(valid.names, ".", fe)
  names(valid.urls) <- valid.filenames
  return(list(df = valid.df, urls = valid.urls))
}

# used in hca, recursively extract projects (limit size to 100, get all projects when the size is greater than 100)
# reference: https://bioconductor.org/packages/release/bioc/html/hca.html
RecurURLRetrieval <- function(url) {
  url.content <- URLRetrieval(url)
  next.url <- url.content$pagination$`next`
  if (!is.null(next.url)) {
    # return(c(url.content, RecurURLRetrieval(next.url)))
    return(data.table::rbindlist(list(url.content$hits, RecurURLRetrieval(next.url)), fill = TRUE))
  } else {
    return(url.content$hits)
  }
}

# used in hca, two-level list, final is vector
HCAPasteCol <- function(df, col) {
  if (col %in% colnames(df)) {
    col.value <- paste0(sapply(
      df[[col]],
      function(x) {
        ifelse(is.null(x), "",
          paste0(x, collapse = "|")
        )
      }
    ), collapse = ", ")
  } else {
    col.value <- ""
  }
  return(col.value)
}

# used in hca, dataframe, check column exists
HCAPasteColdf <- function(df, col = NULL) {
  if (col %in% colnames(df)) {
    return(paste0(df[[col]], collapse = ", "))
  } else {
    return("")
  }
}

# used in hca, filter proojects
HCAAttrFilter <- function(df, attr, attr.value) {
  if (is.null(attr.value)) {
    message("Use all ", attr, " as input!")
    # return row index
    value <- 1:nrow(df)
  } else {
    # lower case
    attr.value <- tolower(attr.value)
    all.values <- df[[attr]]
    # value contains dot
    value.list <- list()
    for (pv in attr.value) {
      pv.vec <- c()
      for (avi in 1:length(all.values)) {
        av <- all.values[avi]
        av.vec <- strsplit(x = strsplit(x = av, split = ", ")[[1]], split = "\\|") %>%
          unlist() %>%
          tolower()
        if (pv %in% av.vec) {
          pv.vec <- c(pv.vec, avi)
        }
      }
      value.list[[pv]] <- pv.vec
    }
    # get invalid value
    invalid.value <- setdiff(attr.value, names(value.list))
    value <- unique(unlist(value.list))
    # print invalid value
    if (length(invalid.value) > 0) {
      message(paste0(invalid.value, collapse = ", "), " are not valid for ", attr)
    }
    # deal with empty value
    if (length(value) == 0) {
      warning("There is no valid value under ", paste(attr.value, collapse = ", "), " in ", attr, ". This filter is invalid!")
      value <- 1:nrow(df)
    }
  }
  return(value)
}

# get filter value for "PanglaoDB", "UCSC", "CELLxGENE", "HCA"
CheckFilter <- function(df, filter, all.filter, database, combine) {
  if (combine) {
    # extract dataframe
    filter.df <- df[, all.filter[filter]]
    # filter all empty row
    filter.df <- dplyr::filter(filter.df, !dplyr::if_all(
      dplyr::everything(),
      function(x) {
        x == "" | is.na(x)
      }
    ))
    # for UCSC, remove parent
    if (database == "UCSC") {
      filter.df <- apply(filter.df, 2, FUN = function(x) {
        gsub("\\|parent$", replacement = "", x = tolower(x))
      }) %>% as.data.frame()
    }
    # split dataframe, a dataset may contain multiple value, eg: multiple organ
    ## sep is ,
    for (col in colnames(filter.df)) {
      filter.df <- tidyr::separate_rows(filter.df, tidyr::all_of(col), sep = ", ")
    }
    ## sep is |
    for (col in colnames(filter.df)) {
      filter.df <- tidyr::separate_rows(filter.df, tidyr::all_of(col), sep = "\\|")
    }
    if (database != "PanglaoDB") {
      filter.df <- apply(filter.df, 2, tolower) %>% as.data.frame()
    }
    # summarise
    filter.df.stat <- filter.df %>%
      dplyr::group_by_all() %>%
      dplyr::summarise(Num = dplyr::n()) %>%
      dplyr::arrange(dplyr::desc(Num))
    return(filter.df.stat)
  } else {
    filter.list <- lapply(filter, function(x) {
      filter.values <- df[[all.filter[x]]]
      if (database == "UCSC") {
        filter.values <- gsub("\\|parent$", replacement = "", x = tolower(filter.values))
      }
      if (database == "PanglaoDB") {
        vf.df <- strsplit(x = unlist(strsplit(x = filter.values, split = ", ")), split = "\\|") %>%
          unlist() %>%
          table() %>%
          as.data.frame()
      } else {
        vf.df <- strsplit(x = unlist(strsplit(x = filter.values, split = ", ")), split = "\\|") %>%
          unlist() %>%
          tolower() %>%
          table() %>%
          as.data.frame()
      }
      colnames(vf.df) <- c("Value", "Num")
      vf.df <- vf.df[order(vf.df$Num, decreasing = TRUE), ]
      vf.df$Key <- x
      rownames(vf.df) <- NULL
      return(vf.df)
    })
    names(filter.list) <- filter
    return(filter.list)
  }
}

# used in hca, extract data information from contributedAnalyses and matrices
HCAExtactData <- function(df) {
  # unlist
  df.vec <- unlist(df)
  # create dataframe
  df.unlist <- data.frame(meta = names(df.vec), value = df.vec)
  rownames(df.unlist) <- NULL
  # data columns
  data.cols <- c(
    "contentDescription", "format", "isIntermediate", "name", "sha256", "size",
    "fileSource", "uuid", "version", "matrixCellCount", "drs_uri", "url"
  )
  data.col.pattern <- paste0(data.cols, collapse = "|")
  type.pattern <- paste0("(.*)\\.(", data.col.pattern, ")([0-9]*)")
  # add col
  df.unlist$type <- gsub(pattern = type.pattern, replacement = "\\2", x = df.unlist$meta)
  df.unlist$num <- gsub(pattern = type.pattern, replacement = "\\3", x = df.unlist$meta)
  df.unlist$meta <- gsub(pattern = type.pattern, replacement = "\\1", x = df.unlist$meta)
  df.unlist$meta <- paste0(df.unlist$meta, ".", df.unlist$num)
  df.final <- tidyr::spread(data = df.unlist[c("meta", "type", "value")], key = "type", value = "value")
  return(df.final)
}

# used in CELLxGENE, Zenodo
LoadRDS2Seurat <- function(out.folder, merge, obs.value.filter = NULL, obs.keys = NULL, include.genes = NULL) {
  rds.files <- list.files(path = out.folder, pattern = "rds$", full.names = TRUE, ignore.case = TRUE)
  if (length(rds.files) > 0) {
    message("There is rds in file.ext and return.seu is TRUE, return SeuratOnject!")
    seu.list <- sapply(X = rds.files, FUN = function(x) {
      tryCatch(
        {
          x.rds <- readRDS(x)
          if (class(x.rds) == "Seurat") {
            if (!is.null(obs.value.filter) || !is.null(obs.keys) || !is.null(include.genes)) {
              x.rds.df <- x.rds@meta.data
              # filter dataset metadata
              ## filter cell's metadata values
              if (!is.null(obs.value.filter)) {
                x.rds.df <- x.rds.df %>% dplyr::filter(eval(rlang::parse_expr(obs.value.filter)))
              }
              ## filter cell's metadata colnames
              if (!is.null(obs.keys)) {
                x.rds.df <- x.rds.df[obs.keys]
              }
              x.rds <- subset(x = x.rds, cells = rownames(x.rds.df), features = include.genes)
            }
            x.rds
          } else {
            message(x, " is not SeuratObject, skip!")
            NULL
          }
        },
        error = function(cond) {
          message("Reading ", x, " error:", cond)
          NULL
        }
      )
    })
    if (isTRUE(merge)) {
      seu.obj <- mergeExperiments(seu.list)
    } else {
      seu.obj <- seu.list
    }
    return(seu.obj)
  } else {
    message("There is no rds file under ", out.folder)
    return(NULL)
  }
}

checkGeneNames <- function(seurat, n = 6) {
  gene_names <- seurat %>% rownames() %>% head(n)
  log_info("Please check if the gene names:{gene_names} are correct.")
  # return(gene_names)
}

checkSeuratAssay <- function(seurat,nrow = 6,ncol=4) {
  # log_info("Extracting first {nrow} rows and {ncol} columns of assay data from Seurat object.")
  log_info("Extracting several rows and columns of assay data from Seurat object.")
  
  data <- GetAssayData(seurat)[1:nrow, 1:ncol] %>% 
    as.data.frame() %>%
    rownames_to_column("gene")
  
  log_info(data)
  
  log_info("Please check if the gene and cell names are correct.")
  
}

checkFileExtension <- function(file_path,  allowed_extensions = c("rds", "rdata", "h5ad", "loom")) {
  # Define allowed file extensions
  # allowed_extensions <- c("rds", "Rdara", "h5ad", "loom")
  compressed_extensions <- c("gz", "tar")
  
  # Extract file extension
  file_extension <- tools::file_ext(file_path)
  
  # Check if the file extension is a compressed file type
  if (tolower(file_extension) %in% compressed_extensions) {
    # Remove compressed extension and check the actual file extension
    base_name <- tools::file_path_sans_ext(file_path)
    base_extension <- tools::file_ext(base_name)
    is_valid <- tolower(base_extension) %in% tolower(allowed_extensions)
  } else {
    # Check if the file extension is within the allowed range
    is_valid <- tolower(file_extension) %in% tolower(allowed_extensions)
  }
  
  return(is_valid)
}


# Find OBJ file from all file paths from supp.type given
findOBJFile = function(supp.type, filePaths) {
  typeExtensions <- list(
    Seurat = c("rds", "rdata"),
    SingleCellExperiment = c("rds", "rdata"),
    Loom = c("loom"),
    AnnData = c("h5ad")
  )
  
  if (!supp.type %in% names(typeExtensions)) {
    log_error("Invalid supp.type provided.")
    return(NULL)
  }
  
  allowedExtensions <- typeExtensions[[supp.type]]
  validPaths <- c()
  
  for (filePath in filePaths) {
    # Extract file extension
    fileExtension <- tools::file_ext(filePath)
    
    # Check if the file extension is a compressed file type
    if (fileExtension %in% c("gz", "tar")) {
      # Remove compressed extension and check the actual file extension
      baseName <- tools::file_path_sans_ext(filePath)
      baseExtension <- tools::file_ext(baseName)
      isValid <- tolower(baseExtension) %in% tolower(allowedExtensions)
    } else {
      # Check if the file extension is within the allowed range
      isValid <- tolower(fileExtension) %in% tolower(allowedExtensions)
    }
    
    if (isValid) {
      validPaths <- c(validPaths, filePath)
    }
  }
  
  return(validPaths)
}

processJson <- function(json_file_path) {
  data <- fromJSON(json_file_path)
  
  source_id <- data$source_id
  file_type <- data$file_type
  has_meta = data$meta_file_pattern
  
  supp_idx <- ifelse(is.null(data$supp.idx) || data$supp.idx == "", 2, data$supp.idx)
  
  list(
    source_id = source_id,
    file_type = file_type,
    has_meta = has_meta,
    supp_idx = supp_idx
  )
}

isCountMTX = function(file){
  # Read the first 10 lines of the file using zcat
  con <- gzfile(file, "r")
  lines <- readLines(con, n = 10)
  close(con)
  
  # Skip the header line starting with "%%MatrixMarket"
  data_lines <- lines[-1]
  
  # Split lines into individual values and convert to numeric
  data_values <- unlist(strsplit(data_lines, "\\s+"))
  data_values <- as.numeric(data_values)
  
  # Check if all values are integers
  all_integer <- all(data_values == floor(data_values))
  
  return(all_integer)
}



# Define function to check if a file is a counts matrix
isCountsMatrix <- function(file_path) {
  # Check the file extension
  file_ext <- tools::file_ext(file_path)
  if (file_ext == c("gz")) {
    # Decompress the file
      gunzip(file_path, overwrite = TRUE)
      file_path <- gsub(pattern = "\\.gz", replacement = "", x = file_path)
  }
  # 
  # if (grepl("tar$", file_path)) {
  #   exdir =  dirname(file_path)
  # 
  #   extractTar(file_path,exdir =exdir )
  #   file_path <- gsub(pattern = "\\.tar", replacement = "", x = file_path)
  # }  
  # Use the first decompressed file
  file_ext <- tools::file_ext(file_path)
  
  if (!file_ext %in% c("xlsx", "xls", "csv", "tsv", "txt")) {
    return(FALSE)
  }
  
  # Read the file based on its extension
  if (file_ext %in% c("xlsx", "xls")) {
    log_info("check {file_path}...")
    df_full <- read_excel(file_path)
  } else if (file_ext == "csv") {
    log_info("check {file_path}...")
    df_full <- fread(file_path, sep = ",")
  } else if (file_ext == "tsv" || file_ext == "txt") {
    log_info("check {file_path}...")
    df_full <- fread(file_path, sep = "\t")
  } else {
    return(FALSE) # Return FALSE for unsupported file formats
  }
  
  if(min(dim(df_full))<1000){
    return(FALSE)
  }
  
  # Check if column names contain cell barcodes
  # cell_barcode_check <- all(grepl("^cell_barcode", colnames(df)))
  df = df_full[1:5,1:5]
  log_info(df)
  # Check if row names are gene symbols
  # gene_symbol_check <- all(grepl("^[A-Za-z0-9_-]+$", df[[1]]))
  # log_info(df)
  
  # Remove the gene symbol column
  df_full <- df_full[, -1, with = FALSE]
  
  # Check if all values are integers
  integer_check <- all(sapply(df_full, is.integer))
  
  # Return the result of all checks
  return(integer_check)
}

# Define a function to check if a directory is empty and delete files if not
clean_directory_if_not_empty <- function(directory) {
  # Get all files and directories within the directory
  files <- list.files(directory,  full.names = TRUE, recursive = TRUE)
  
  # Retain only files, exclude directories
  file_paths <- files[file.info(files)$isdir == FALSE]
  
  dir_paths <- list.dirs(directory, recursive = TRUE, full.names = TRUE)

  
  # If the directory is not empty, delete the files and directories
  if (length(file_paths) > 0) {
    file.remove(file_paths)
    message(str_glue("Files in directory {directory} have been deleted."))
  }
  
  # If there are directories, delete them
  if (length(dir_paths) > 0) {
    # Sort directories by depth (deeper directories first) to avoid errors
    dir_paths <- dir_paths[order(nchar(dir_paths), decreasing = TRUE)]
    sapply(dir_paths, unlink, recursive = TRUE, force = TRUE)
    message(str_glue("Subdirectories in directory {directory} have been deleted."))
  }
  
  # Final message if the directory is empty
  if (length(file_paths) == 0 && length(dir_paths) == 0) {
    message(str_glue("Directory {directory} is empty."))
  }
}


extractTar <- function(file_path,exdir = NULL) {
  if(grepl(".tar.gz$",file_path)){
    gunzip(file_path,overwrite  = T)
    file_path <- gsub(pattern = "\\.gz", replacement = "", x = file_path)
  }
  # Create the extraction directory path
  if (is.null(exdir)) {
    file_dir <- dirname(file_path)
    exdir <- file.path(file_dir, digest::digest(file_path, algo = "md5"))
  }else{
    exdir = exdir
  }
  # Ensure the extraction directory exists
  if (!dir.exists(exdir)) {
    dir.create(path = exdir, recursive = TRUE)
  }
  
  # Try to extract the .tar file to the extraction directory
  tryCatch({
    utils::untar(file_path, exdir = exdir)
    log_info(str_glue("File {file_path} has been extracted to {exdir}."))
    file.remove(file_path)
  }, warning = function(w) {
    log_warn(str_glue("Warning during extraction of {file_path}: {w$message}"))
    warning(w)
  }, error = function(e) {
    log_error(str_glue("Error during extraction of {file_path}: {e$message}"))
    stop(e)
  })
}

findCommonPart <- function(paths) {
  # Helper function to find the common prefix of two strings
  common_prefix <- function(x, y) {
    str_common <- ""
    for(i in seq_len(min(nchar(x), nchar(y)))) {
      if(substr(x, i, i) == substr(y, i, i)) {
        str_common <- paste0(str_common, substr(x, i, i))
      } else {
        break
      }
    }
    return(str_common)
  }
    common_part <- Reduce(common_prefix, paths)
    
    # Replace paths with the common part
    new_paths <-rep(common_part,length(paths))
    names(new_paths) = names(paths)
    
    # Return the new paths
    return(new_paths)
}


findLargestMatrix <- function(matrix_list) {
  sizes <- sapply(matrix_list, function(x) prod(dim(x)))
  largest_matrix <- matrix_list[[which.max(sizes)]]
  return(largest_matrix)
}



checkBarcodes <- function(file_path) {
  data <- fread(file_path)
  
  num_cols <- ncol(data)
  
  if (num_cols != 1) {
    log_warn("The file has multiple columns. Identifying the most likely barcode column.")

    col_lengths <- sapply(data, function(col) sum(nchar(as.character(col))))
    barcode_col <- which.max(col_lengths)
    
    barcodes <- data %>% select(barcode_col)
    log_info(head(barcodes))
    write.table(barcodes, file = gzfile(file_path), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
    }
  
  return()
}
