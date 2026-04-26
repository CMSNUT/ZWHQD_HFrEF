# Clean function: ensure CID contains only digits
clean_cids <- function(cid_vec) {
  cid_vec <- as.character(cid_vec)
  # Remove all spaces, tabs, newlines
  cid_vec <- gsub("[[:space:]]", "", cid_vec)
  # Keep only digits (if letters or other characters, set to NA)
  cid_vec <- ifelse(grepl("^[0-9]+$", cid_vec), cid_vec, NA)
  return(cid_vec[!is.na(cid_vec)])
}

# Define fetch_batch function, always return data.frame (even on failure)
fetch_batch <- function(properties, batch_cids, retries = 10, delay = 0.5) {
  # 1. Clean input
  # Clean current batch CIDs
  batch_cids <- clean_cids(batch_cids)
  
  if (length(batch_cids) == 0) return(data.frame())
  
  # 2. Build request
  all_cids_str <- paste(batch_cids, collapse = ",")
  api_url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                    all_cids_str, "/property/", paste(properties, collapse = ","), "/JSON")
  
  # 3. Retry loop to handle rate limiting
  for (attempt in seq_len(retries)) {
    result <- tryCatch({
      # Add delay before each request
      Sys.sleep(delay)
      res <- httr::GET(api_url)
      
      # Handle different HTTP status codes
      if (httr::status_code(res) == 200) {
        data <- jsonlite::fromJSON(httr::content(res, "text", encoding = "UTF-8"))
        df <- data$PropertyTable$Properties
        if (!is.null(df) && nrow(df) > 0) {
          df$CID <- as.character(df$CID)
          return(df)
        } else {
          return(data.frame())
        }
      } else if (httr::status_code(res) == 503) {
        # Server overload, wait longer and retry
        warning(sprintf("Batch %s encountered 503 error, retry %d...", all_cids_str, attempt))
        Sys.sleep(delay * 2^attempt) # exponential backoff
        next
      } else {
        warning(sprintf("Batch failed, status code: %s", httr::status_code(res)))
        return(data.frame())
      }
    }, error = function(e) {
      warning(sprintf("Batch request error: %s", e$message))
      return(NULL)
    })
    
    # If data successfully retrieved, return immediately
    if (!is.null(result) && nrow(result) > 0) {
      return(result)
    }
  }
  
  # All retries failed, return empty data frame
  warning(sprintf("Batch %s ultimately failed.", all_cids_str))
  return(data.frame())
}

# read_compound_target
read_compound_target <- function(folder) {
  # Check if folder exists
  if (!dir.exists(folder)) {
    stop("Folder does not exist: ", folder)
  }
  
  # Get all csv and xlsx files
  csv_files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE)
  xlsx_files <- list.files(folder, pattern = "\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
  files <- c(csv_files, xlsx_files)
  
  if (length(files) == 0) {
    warning("No CSV or XLSX files found in the folder.")
    return(data.frame())
  }
  
  # Load required packages (error if not installed)
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required to read .xlsx files. Please install it.")
  }
  
  # Define function to read a single file
  read_one_file <- function(file_path) {
    # Get filename (without extension) as CID
    cid <- tools::file_path_sans_ext(basename(file_path))
    
    # Choose reading method based on extension
    ext <- tolower(tools::file_ext(file_path))
    if (ext == "csv") {
      df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (ext == "xlsx") {
      df <- readxl::read_excel(file_path, .name_repair = "minimal")  # keep original column names
      df <- as.data.frame(df, stringsAsFactors = FALSE)
    } else {
      stop("Unsupported file type: ", ext)
    }
    
    # Add CID column as first column
    df <- cbind(CID = cid, df)
    return(df)
  }
  
  # Read all files and combine
  all_dfs <- lapply(files, read_one_file)
  combined <- do.call(rbind, lapply(all_dfs, function(df) {
    # Using dplyr::bind_rows automatically matches column names and fills NAs, but to avoid dependencies, use rbind.fill or base method
    # Here use data.table::rbindlist for efficiency, but to minimize dependencies, implement simple column-name merge
    df
  }))
  
  # Use rbind.fill method for inconsistent columns (from plyr or simple implementation)
  # To avoid extra dependencies, use dplyr::bind_rows if available, otherwise base method
  if (requireNamespace("dplyr", quietly = TRUE)) {
    combined <- dplyr::bind_rows(all_dfs)
  } else {
    # Manual merge: get union of all column names, fill missing columns with NA
    all_cols <- unique(unlist(lapply(all_dfs, names)))
    combined <- do.call(rbind, lapply(all_dfs, function(df) {
      missing <- setdiff(all_cols, names(df))
      for (col in missing) df[[col]] <- NA
      df[, all_cols, drop = FALSE]
    }))
  }
  
  return(combined)
}

get_symbol_by_uniprot <- function(uniprot_ids, unique_by = "UNIPROT") {
  # Keep only Swiss-Prot formatted input IDs
  swissprot_pattern <- "^[OPQ][0-9][A-Z0-9]{3}[0-9]($|[0-9])"
  uniprot_ids <- uniprot_ids[grepl(swissprot_pattern, uniprot_ids)]
  
  uniprot_clean <- sub("\\.[0-9]+$", "", uniprot_ids)
  uniprot_clean <- unique(uniprot_clean)
  
  # Use multiVals = "first" to return only the first match per key
  result <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = uniprot_clean,
                                  keytype = "UNIPROT",
                                  columns = c("SYMBOL", "UNIPROT"),
                                  multiVals = "first") %>%   # key change
    dplyr::select(Symbol = SYMBOL, Uniprot.ID = UNIPROT) %>%
    dplyr::mutate(dplyr::across(where(is.character), trimws))
  
  # If first_only = TRUE, multiVals = "first" already handles it, no need to repeat.
  # However, keep deduplication by unique_by if needed:
  if (unique_by == "SYMBOL") {
    result <- result %>% dplyr::distinct(Symbol, .keep_all = TRUE)
  } else if (unique_by == "UNIPROT") {
    result <- result %>% dplyr::distinct(Uniprot.ID, .keep_all = TRUE)
  } else {
    stop("unique_by must be 'SYMBOL' or 'UNIPROT'")
  }
  
  return(result)
}

get_uniprot_by_symbol <- function(symbols, unique_by = "SYMBOL") {
  symbols <- unique(symbols)
  
  # Use multiVals = "first" to return only the first UNIPROT per Symbol
  all_mappings <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = symbols,
                                        columns = c("SYMBOL", "UNIPROT"),
                                        keytype = "SYMBOL",
                                        multiVals = "first")
  
  # Keep only Swiss-Prot format
  swissprot_pattern <- "^[OPQ][0-9][A-Z0-9]{3}[0-9]($|[0-9])"
  all_mappings_swiss <- all_mappings %>%
    dplyr::filter(grepl(swissprot_pattern, UNIPROT))
  
  # Deduplicate by unique_by
  if (unique_by == "SYMBOL") {
    result <- all_mappings_swiss %>%
      dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
      dplyr::select(Symbol = SYMBOL, Uniprot.ID = UNIPROT)
  } else if (unique_by == "UNIPROT") {
    result <- all_mappings_swiss %>%
      dplyr::distinct(UNIPROT, .keep_all = TRUE) %>%
      dplyr::select(Symbol = SYMBOL, Uniprot.ID = UNIPROT)
  } else {
    stop("unique_by must be 'SYMBOL' or 'UNIPROT'")
  }
  
  return(result)
}
