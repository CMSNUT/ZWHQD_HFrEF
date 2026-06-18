#==============================================================================
# Clean CID: keep only numeric CIDs, remove invalid/empty values
#==============================================================================
clean_cids <- function(cid_vec) {
  cid_vec <- as.character(cid_vec)
  cid_vec <- gsub("[[:space:]]", "", cid_vec)  # Remove all whitespace
  cid_vec <- ifelse(grepl("^[0-9]+$", cid_vec), cid_vec, NA)  # Keep digits only
  return(cid_vec[!is.na(cid_vec)])
}

#==============================================================================
# Batch query PubChem API with retry & exponential backoff
#==============================================================================
fetch_batch <- function(properties, batch_cids, retries = 10, delay = 0.5) {
  # Clean input CIDs
  batch_cids <- clean_cids(batch_cids)
  if (length(batch_cids) == 0) return(data.frame())
  
  # Build API URL
  cid_str <- paste(batch_cids, collapse = ",")
  prop_str <- paste(properties, collapse = ",")
  api_url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                    cid_str, "/property/", prop_str, "/JSON")
  
  # Retry loop
  for (attempt in seq_len(retries)) {
    result <- tryCatch({
      Sys.sleep(delay)
      res <- httr::GET(api_url)
      
      # Success
      if (httr::status_code(res) == 200) {
        data <- jsonlite::fromJSON(httr::content(res, "text", encoding = "UTF-8"))
        df <- data$PropertyTable$Properties
        if (!is.null(df) && nrow(df) > 0) {
          df$CID <- as.character(df$CID)
          return(df)
        } else {
          return(data.frame())
        }
      }
      
      # Server overload: exponential backoff
      if (httr::status_code(res) == 503) {
        warning(sprintf("503 | CID batch: %s | retry %s", cid_str, attempt))
        Sys.sleep(delay * 2^attempt)
        next
      }
      
      # Other HTTP errors
      warning(sprintf("Request failed | code: %s", httr::status_code(res)))
      return(data.frame())
      
    }, error = function(e) {
      warning(sprintf("Error: %s", e$message))
      return(NULL)
    })
    
    if (!is.null(result) && nrow(result) > 0) return(result)
  }
  
  warning(sprintf("All retries failed | CID batch: %s", cid_str))
  return(data.frame())
}


#==============================================================================
# download SDF 
#==============================================================================

download_all_sdf <- function(cid_list, output_dir = "sdf_files") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  results <- data.frame()
  total <- length(cid_list)
  
  for (i in seq_along(cid_list)) {
    cid <- cid_list[i]
    cat(sprintf("[%d/%d] Downloading CID %s ... ", i, total, cid))
    
    url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid, "/SDF")
    
    # Try 3D
    res <- httr::GET(url, query = list(record_type = "3d"))
    
    if (status_code(res) == 200) {
      out_file <- file.path(output_dir, paste0(cid, ".sdf"))
      writeBin(content(res, "raw"), out_file)
      cat("OK (3D)\n")
      results <- rbind(results, data.frame(CID = cid, status = "success", type = "3d"))
      next
    }
    
    # Try 2D
    res <- httr::GET(url, query = list(record_type = "2d"))
    if (status_code(res) == 200) {
      out_file <- file.path(output_dir, paste0(cid, ".sdf"))
      writeBin(content(res, "raw"), out_file)
      cat("OK (2D)\n")
      results <- rbind(results, data.frame(CID = cid, status = "success", type = "2d"))
      next
    }
    
    cat("FAILED\n")
    results <- rbind(results, data.frame(CID = cid, status = "failed", type = "none"))
    
    Sys.sleep(0.2)
  }
  
  # Summary
  cat("\n========== Summary ==========\n")
  cat("Total:", nrow(results), "\n")
  cat("Success (3D):", sum(results$type == "3d"), "\n")
  cat("Success (2D):", sum(results$type == "2d"), "\n")
  cat("Failed:", sum(results$status == "failed"), "\n")
  
  write.csv(results, file.path(output_dir, "download_report.csv"), row.names = FALSE)
  
  return(results)
}



#==============================================================================
# Read compound-target files (CSV/XLSX) from folder, add CID from filename
#==============================================================================
read_compound_target <- function(folder, skip_row = 0) {
  if (!dir.exists(folder)) stop("Folder not found: ", folder)
  
  # List all files
  csv_files  <- list.files(folder, "\\.csv$", full.names = TRUE, ignore.case = TRUE)
  xlsx_files <- list.files(folder, "\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
  files <- c(csv_files, xlsx_files)

  if (length(files) == 0) {
    warning("No CSV/XLSX files found.")
    return(data.frame())
  }
  
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required for XLSX files.")
  }
  
  # Read single file
  read_one_file <- function(file_path, skip_row) {
    cid <- tools::file_path_sans_ext(basename(file_path))
    ext <- tolower(tools::file_ext(file_path))
    
    if (ext == "csv") {
      df <- read.csv(file_path, 
                     stringsAsFactors = FALSE, 
                     check.names = FALSE, 
                     skip = skip_row,
                     row.names = NULL)
    } else if (ext == "xlsx") {
      df <- as.data.frame(
        readxl::read_excel(file_path, 
                           .name_repair = "minimal", 
                           skip = skip_row),
        stringsAsFactors = FALSE
      )
    } else {
      stop("Unsupported file type: ", ext)
    }
    
    df <- cbind(CID = cid, df)
    return(df)
  }
  
  # Read and merge
  all_dfs <- lapply(files, function(f) read_one_file(f, skip_row))
  
  # Safe row-bind (auto-align columns)
  if (requireNamespace("dplyr", quietly = TRUE)) {
    combined <- dplyr::bind_rows(all_dfs)
  } else {
    all_cols <- unique(unlist(lapply(all_dfs, names)))
    combined <- do.call(rbind, lapply(all_dfs, function(df) {
      missing <- setdiff(all_cols, names(df))
      for (col in missing) df[[col]] <- NA
      df[, all_cols, drop = FALSE]
    }))
  }
  
  return(combined)
}


read_compound_target2 <- function(folder, skip_row = 0) {
  if (!dir.exists(folder)) stop("Folder not found: ", folder)
  
  # List all files
  csv_files  <- list.files(folder, "\\.csv$", full.names = TRUE, ignore.case = TRUE)
  xlsx_files <- list.files(folder, "\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
  files <- c(csv_files, xlsx_files)
  
  if (length(files) == 0) {
    warning("No CSV/XLSX files found.")
    return(data.frame())
  }
  
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required for XLSX files.")
  }
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required for CSV files.")
  }
  
  # Read single file
  read_one_file <- function(file_path, skip_row) {
    cid <- tools::file_path_sans_ext(basename(file_path))
    ext <- tolower(tools::file_ext(file_path))
    
    if (ext == "csv") {
      df <- readr::read_csv(file_path, 
                            skip = skip_row,
                            show_col_types = FALSE) 
      df <- as.data.frame(df, stringsAsFactors = FALSE)
    } else if (ext == "xlsx") {
      df <- as.data.frame(
        readxl::read_excel(file_path, 
                           .name_repair = "minimal", 
                           skip = skip_row),
        stringsAsFactors = FALSE
      )
    } else {
      stop("Unsupported file type: ", ext)
    }
    
    df <- cbind(CID = cid, df)
    return(df)
  }
  
  # Read and merge
  all_dfs <- lapply(files, function(f) read_one_file(f, skip_row))
  
  # Safe row-bind (auto-align columns)
  if (requireNamespace("dplyr", quietly = TRUE)) {
    combined <- dplyr::bind_rows(all_dfs)
  } else {
    all_cols <- unique(unlist(lapply(all_dfs, names)))
    combined <- do.call(rbind, lapply(all_dfs, function(df) {
      missing <- setdiff(all_cols, names(df))
      for (col in missing) df[[col]] <- NA
      df[, all_cols, drop = FALSE]
    }))
  }
  
  return(combined)
}


chembl_to_uniprot <- function(chembl_id) {
  url <- paste0("https://www.ebi.ac.uk/chembl/api/data/target/", chembl_id, ".json")
  
  response <- GET(url)
  
  if (http_error(response)) {
    warning(paste("API request error. Status Code:", status_code(response)))
    return(NA)
  }
  
  # JSON
  data <- fromJSON(content(response, "text", encoding = "UTF-8"))
  
  # target_components
  components <- data$target_components
  
  if (is.null(components) || nrow(components) == 0) {
    warning(paste(chembl_id, "没有关联的 UniProt 条目"))
    return(NA)
  }
  
  # get UniProt accession
  uniprot_ids <- unique(components$accession)
  uniprot_ids <- uniprot_ids[!is.na(uniprot_ids)]
  
  if (length(uniprot_ids) == 0) {
    return(NA)
  } else if (length(uniprot_ids) == 1) {
    return(uniprot_ids[1])
  } else {
    return(uniprot_ids)
  }
}


#==============================================================================
# Map Uniprot ID → Gene Symbol (keep only Swiss-Prot IDs)
#==============================================================================
get_symbol_by_uniprot <- function(uniprot_ids, unique_by = "UNIPROT") {
  # Swiss-Prot pattern filter
  swiss_pattern <- "^[OPQ][0-9][A-Z0-9]{3}[0-9]($|[0-9])"
  uniprot_ids <- uniprot_ids[grepl(swiss_pattern, uniprot_ids)]
  uniprot_clean <- sub("\\.[0-9]+$", "", uniprot_ids)
  uniprot_clean <- unique(uniprot_clean)
  
  # Query with 1:1 mapping
  result <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = uniprot_clean,
    keytype = "UNIPROT",
    columns = c("SYMBOL", "UNIPROT"),
    multiVals = "first"
  ) %>%
    dplyr::select(Symbol = SYMBOL, Uniprot.ID = UNIPROT) %>%
    dplyr::mutate(dplyr::across(where(is.character), trimws))
  
  # Deduplication
  if (unique_by == "SYMBOL") {
    result <- dplyr::distinct(result, Symbol, .keep_all = TRUE)
  } else if (unique_by == "UNIPROT") {
    result <- dplyr::distinct(result, Uniprot.ID, .keep_all = TRUE)
  } else {
    stop("unique_by must be 'SYMBOL' or 'UNIPROT'")
  }
  
  return(result)
}

#==============================================================================
# Map Gene Symbol → Uniprot ID (Swiss-Prot only)
#==============================================================================
get_uniprot_by_symbol <- function(symbols, unique_by = "SYMBOL") {
  symbols <- unique(symbols)
  swiss_pattern <- "^[OPQ][0-9][A-Z0-9]{3}[0-9]($|[0-9])"
  
  # 1:1 mapping
  all_mappings <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "UNIPROT"),
    multiVals = "first"
  )
  
  # Keep valid Swiss-Prot entries
  all_mappings_swiss <- dplyr::filter(all_mappings, grepl(swiss_pattern, UNIPROT))
  
  # Deduplication
  if (unique_by == "SYMBOL") {
    result <- dplyr::distinct(all_mappings_swiss, SYMBOL, .keep_all = TRUE)
  } else if (unique_by == "UNIPROT") {
    result <- dplyr::distinct(all_mappings_swiss, UNIPROT, .keep_all = TRUE)
  } else {
    stop("unique_by must be 'SYMBOL' or 'UNIPROT'")
  }
  
  result <- dplyr::select(result, Symbol = SYMBOL, Uniprot.ID = UNIPROT)
  return(result)
}


