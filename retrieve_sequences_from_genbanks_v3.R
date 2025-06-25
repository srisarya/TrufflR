# =============================================================================
# GENE SEQUENCE EXTRACTION FUNCTIONS
# =============================================================================
# This module contains functions for extracting gene sequences from NCBI
# based on taxonomic IDs and gene synonyms. The main workflow:
# 1. Read taxonomic IDs from file
# 2. Search NCBI for sequences matching genes of interest
# 3. Download and parse GenBank/FASTA records
# 4. Extract specific gene regions
# 5. Generate organized output files and summaries

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Extract accession number from FASTA header
#' 
#' @param header_line Character string containing FASTA header
#' @return Character string with accession number
extract_accession <- function(header_line) {
  # Remove the leading ">" and extract accession from various header formats
  header_clean <- gsub("^>", "", header_line)
  
  # Try different patterns to extract accession
  if (grepl("^[A-Z]{1,2}[0-9]{5,8}(\\.[0-9]+)?", header_clean)) {
    # Standard GenBank format (e.g., "AB123456.1")
    accession <- regmatches(header_clean, regexpr("^[A-Z]{1,2}[0-9]{5,8}(\\.[0-9]+)?", header_clean))
  } else {
    # Fallback: take first word before space or pipe
    accession <- strsplit(header_clean, "[ |]")[[1]][1]
  }
  
  return(accession)
}

#' Determine sequence type from title
#' 
#' @param title Character string containing sequence title
#' @return Character string indicating sequence type
determine_sequence_type <- function(title) {
  title_lower <- tolower(title)
  # Define sequence type patterns
  if (grepl("complete genome|whole genome", title_lower)) {
    return("complete_genome")
  } else if (grepl("chromosome|chr", title_lower)) {
    return("chromosome")
  } else if (grepl("plasmid", title_lower)) {
    return("plasmid")
  } else if (grepl("mitochondr|chloroplast", title_lower)) {
    return("organellar")
  } else if (grepl("scaffold|contig", title_lower)) {
    return("assembly")
  } else if (grepl("gene|cds|mrna", title_lower)) {
    return("gene_region")
  } else {
    return("other")
  }
}

#' Read and validate taxonomic IDs from file
#' 
#' @param taxid_file Path to file containing taxonomic IDs (one per line)
#' @return Character vector of valid taxonomic IDs
#' @details Reads taxid file, removes whitespace and empty lines
read_taxids <- function(taxid_file) {
  cat("Reading taxonomic IDs from:", taxid_file, "\n")
  
  # Check if file exists
  if (!file.exists(taxid_file)) {
    stop("Taxid file not found: ", taxid_file)
  }
  
  # Read and clean taxids
  taxids <- readLines(taxid_file)
  cat("Raw lines read:", length(taxids), "\n")
  
  taxids <- trimws(taxids)  # Remove whitespace
  taxids <- taxids[taxids != ""]  # Remove empty lines
  taxids <- taxids[!is.na(taxids)]  # Remove NA values
  
  cat("Valid taxonomic IDs found:", length(taxids), "\n")
  cat("Taxonomic IDs:", paste(taxids, collapse = ", "), "\n")
  
  if (length(taxids) == 0) {
    stop("No valid taxonomic IDs found in file")
  }
  
  return(taxids)
}

#' Get taxonomic name for a given taxonomic ID
#' 
#' @param taxid Character string containing taxonomic ID
#' @return List containing taxon_name and taxon_name_clean
get_taxon_info <- function(taxid) {
  cat("Retrieving taxonomic information for taxid:", taxid, "\n")
  
  taxon_name <- "Unknown"
  taxon_name_clean <- paste0("taxid_", taxid)  # Fallback name
  
  tryCatch({
    # Fetch taxonomy summary from NCBI
    tax_summary <- rentrez::entrez_summary(db = "taxonomy", id = taxid)
    
    if (!is.null(tax_summary$scientificname)) {
      taxon_name <- tax_summary$scientificname
      
      # Clean taxon name for file naming (remove spaces and special characters)
      taxon_name_clean <- gsub("[^A-Za-z0-9_]", "_", taxon_name)
      taxon_name_clean <- gsub("_+", "_", taxon_name_clean)  # Replace multiple underscores
      taxon_name_clean <- gsub("^_|_$", "", taxon_name_clean)  # Remove leading/trailing underscores
      
      cat("Retrieved taxon name:", taxon_name, "\n")
    }
  }, error = function(e) {
    cat("WARNING: Could not retrieve taxon name for taxid", taxid, ":", e$message, "\n")
  })
  
  return(list(
    taxon_name = taxon_name,
    taxon_name_clean = taxon_name_clean
  ))
}

#' Search NCBI for sequences matching taxonomic ID and gene synonyms
#' 
#' @param taxid Character string containing taxonomic ID
#' @param gene_synonyms Character vector of gene names/synonyms to search for
#' @param retmax Integer maximum number of sequences to retrieve
#' @return List containing search results and query information
search_ncbi_sequences <- function(taxid, gene_synonyms, retmax = 5) {
  cat("Searching NCBI for taxid:", taxid, "\n")
  cat("Gene synonyms:", paste(gene_synonyms, collapse = ", "), "\n")
  
  # Build search query
  gene_query_part <- paste(gene_synonyms, collapse = " OR ")
  full_query <- paste0("txid", taxid, "[Organism] AND (", gene_query_part, ")")
  
  cat("Full NCBI query:", full_query, "\n")
  
  # Search NCBI nuccore database
  search_result <- rentrez::entrez_search(
    db = "nuccore",
    term = full_query,
    retmax = retmax
  )
  
  ids <- search_result$ids
  cat("NCBI search returned", length(ids), "sequence IDs\n")
  
  if (length(ids) > 0) {
    cat("Sequence IDs found:", paste(ids, collapse = ", "), "\n")
  }
  
  return(list(
    ids = ids,
    query = full_query,
    count = length(ids)
  ))
}

#' Create directory structure for output files
#' 
#' @param output_base_dir Base directory for all output files
#' @return List of created directory paths
create_output_directories <- function(output_base_dir) {
  cat("Creating output directory structure in:", output_base_dir, "\n")
  
  # Main output folder
  output_folder <- file.path(output_base_dir)
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Raw files subfolder
  raw_files_folder <- file.path(output_folder, "raw_files")
  if (!dir.exists(raw_files_folder)) {
  }
  
  # Extracted sequences folder (same as main for now)
  extracted_folder <- file.path(output_folder)
  if (!dir.exists(extracted_folder)) {
    dir.create(extracted_folder)
  }
  
  return(list(
    output_folder = output_folder,
    raw_files_folder = raw_files_folder,
    extracted_folder = extracted_folder
  ))
}

#' Download and process sequence records from NCBI
#' 
#' @param ids Character vector of NCBI sequence IDs
#' @param raw_files_folder Directory to store raw GenBank and FASTA files
#' @return List containing processed sequence data and metadata
download_and_process_sequences <- function(ids, raw_files_folder) {
  cat("Processing", length(ids), "sequence records\n")
  
  # Initialize storage structures
  gb_str_dfs <- list()          # GenBank feature dataframes
  gb_seq_dfs <- list()          # FASTA sequence objects  
  sequence_info <- list()       # Sequence metadata
  id_to_accession <- list()     # NCBI ID to accession mapping
  
  # Process each sequence ID
  for (i in seq_along(ids)) {
    id <- ids[i]
    cat("Processing sequence", i, "of", length(ids), "- ID:", id, "\n")
    
    tryCatch({
      # Fetch FASTA record to get accession and title
      cat("Fetching FASTA record for ID:", id, "\n")
      fasta_records <- rentrez::entrez_fetch(
        db = "nuccore",
        id = id,
        rettype = "fasta",
        retmode = "text"
      )
      
      # Extract information from FASTA header
      fasta_lines <- strsplit(fasta_records, "\n")[[1]]
      header_line <- fasta_lines[1]
      sequence_title <- gsub("^>", "", header_line)
      accession <- extract_accession(header_line)
      
      cat("Extracted accession:", accession, "\n")
      cat("Sequence title:", substr(sequence_title, 1, 100), "...\n")  # Truncated for readability
      
      # Store ID to accession mapping
      id_to_accession[[id]] <- accession
      
      # Fetch GenBank record
      cat("Fetching GenBank record for ID:", id, "\n")
      genbank_records <- rentrez::entrez_fetch(
        db = "nuccore",
        id = id,
        rettype = "gb",
        retmode = "text"
      )
      
      # Save files using accession for consistent naming
      filename_gb <- file.path(raw_files_folder, paste0(accession, ".gb"))
      filename_fa <- file.path(raw_files_folder, paste0(accession, ".fa"))
      
      writeLines(genbank_records, filename_gb)
      writeLines(fasta_records, filename_fa)
      
      cat("Saved GenBank file:", filename_gb, "\n")
      cat("Saved FASTA file:", filename_fa, "\n")
      
      # Store sequence metadata
      sequence_info[[id]] <- list(
        id = id,
        accession = accession,
        title = sequence_title,
        type = determine_sequence_type(sequence_title)
      )
      
      # Parse GenBank file for features
      tryCatch({
        cat("Parsing GenBank features for ID:", id, "\n")
        gb_df <- geneviewer::read_gbk(filename_gb) %>% 
          geneviewer::gbk_features_to_df()
        gb_str_dfs[[id]] <- gb_df
        cat("Parsed", nrow(gb_df), "features from GenBank file\n")
      }, error = function(e) {
        cat("WARNING: Error parsing GenBank features for ID", id, ":", e$message, "\n")
      })
      
      # Load FASTA sequence
      tryCatch({
        cat("Loading FASTA sequence for ID:", id, "\n")
        gb_fa <- Biostrings::readDNAStringSet(filename_fa)
        gb_seq_dfs[[id]] <- gb_fa
        cat("Loaded sequence of length:", length(gb_fa[[1]]), "bp\n")
      }, error = function(e) {
        cat("WARNING: Error loading FASTA sequence for ID", id, ":", e$message, "\n")
      })
      
    }, error = function(e) {
      cat("ERROR: Failed to process sequence ID", id, ":", e$message, "\n")
    })
  }
  
  cat("Successfully processed", length(gb_str_dfs), "GenBank files and", 
      length(gb_seq_dfs), "FASTA sequences\n")
  
  return(list(
    gb_str_dfs = gb_str_dfs,
    gb_seq_dfs = gb_seq_dfs,
    sequence_info = sequence_info,
    id_to_accession = id_to_accession
  ))
}

#' Filter GenBank features for target genes and feature types
#' 
#' @param gb_str_dfs List of GenBank feature dataframes
#' @param gene_synonyms Character vector of gene names to search for
#' @param feature_type Character string specifying feature type filter ("all" for no filter)
#' @return List of filtered coordinate dataframes
filter_gene_coordinates <- function(gb_str_dfs, gene_synonyms, feature_type = "all") {
  cat("Filtering gene coordinates for", length(gene_synonyms), "gene synonyms\n")
  cat("Feature type filter:", feature_type, "\n")
  
  # Prepare gene patterns (remove NCBI field tags and convert to lowercase)
  gene_patterns <- gsub("\\[.*?\\]", "", gene_synonyms)
  gene_patterns <- tolower(gene_patterns)
  
  cat("Gene search patterns:", paste(gene_patterns, collapse = ", "), "\n")
  
  # Filter each GenBank dataframe
  gene_coords <- lapply(names(gb_str_dfs), function(id) {
    df <- gb_str_dfs[[id]]
    
    if (!is.null(df) && nrow(df) > 0 && "gene" %in% colnames(df)) {
      cat("Processing", nrow(df), "features for sequence ID:", id, "\n")
      
      # Filter by gene name using pattern matching
      gene_matches <- grepl(paste(gene_patterns, collapse = "|"), 
                            tolower(df$gene), ignore.case = TRUE)
      filtered_df <- df[gene_matches, ]
      
      cat("Found", sum(gene_matches), "gene matches for ID:", id, "\n")
      
      # Filter by feature type if specified
      if (feature_type != "all" && "type" %in% colnames(filtered_df) && nrow(filtered_df) > 0) {
        type_matches <- tolower(filtered_df$type) == tolower(feature_type)
        filtered_df <- filtered_df[type_matches, ]
        cat("After feature type filter:", nrow(filtered_df), "features remain\n")
      }
      
      # Select relevant columns
      available_cols <- intersect(c("gene", "strand", "start", "end", "type"), colnames(filtered_df))
      if (length(available_cols) > 0 && nrow(filtered_df) > 0) {
        result_df <- filtered_df[, available_cols, drop = FALSE]
        cat("Returning", nrow(result_df), "filtered coordinates with columns:", 
            paste(available_cols, collapse = ", "), "\n")
        return(result_df)
      }
    }
    
    cat("No matching features found for sequence ID:", id, "\n")
    return(data.frame())
  })
  
  # Name the list elements and remove empty dataframes
  names(gene_coords) <- names(gb_str_dfs)
  gene_coords <- gene_coords[sapply(gene_coords, nrow) > 0]
  
  cat("Gene coordinates found for", length(gene_coords), "sequences\n")
  
  return(gene_coords)
}

#' Extract gene sequences and save to FASTA files
#' 
#' @param gene_coords List of coordinate dataframes
#' @param gb_seq_dfs List of DNA sequence objects
#' @param sequence_info List of sequence metadata
#' @param id_to_accession Mapping from NCBI IDs to accessions
#' @param extracted_folder Output directory for extracted sequences
#' @param taxon_name_clean Clean taxon name for file naming
#' @param taxid Taxonomic ID
#' @return Number of FASTA files created
extract_and_save_gene_sequences <- function(gene_coords, gb_seq_dfs, sequence_info, 
                                            id_to_accession, extracted_folder, 
                                            taxon_name_clean, taxid) {
  cat("Extracting gene sequences for", length(gene_coords), "sequences\n")
  
  files_created <- 0
  
  # Process each sequence with gene coordinates
  for (id in names(gene_coords)) {
    if (id %in% names(gb_seq_dfs)) {
      coords <- gene_coords[[id]]
      dna_seq <- gb_seq_dfs[[id]]
      accession <- id_to_accession[[id]]
      
      cat("Processing", nrow(coords), "gene regions for accession:", accession, "\n")
      
      # Initialize storage for this sequence's extracted genes
      id_extracted_sequences <- list()
      
      # Extract each gene region
      for (j in 1:nrow(coords)) {
        gene_name <- coords$gene[j]
        start_pos <- coords$start[j]
        end_pos <- coords$end[j]
        strand <- if("strand" %in% colnames(coords)) coords$strand[j] else "+"
        feature_info <- if("type" %in% colnames(coords)) coords$type[j] else "unknown"
        
        cat("Extracting", gene_name, "from positions", start_pos, "-", end_pos, 
            "on strand", strand, "\n")
        
        tryCatch({
          if (length(dna_seq) > 0) {
            full_seq <- dna_seq[[1]]
            
            # Validate coordinates
            if (start_pos > 0 && end_pos <= length(full_seq) && start_pos <= end_pos) {
              gene_subseq <- Biostrings::subseq(full_seq, start = start_pos, end = end_pos)
              
              # Handle reverse complement for minus strand
              if (strand == "-" || strand == -1) {
                gene_subseq <- Biostrings::reverseComplement(gene_subseq)
                cat("Applied reverse complement for minus strand\n")
              }
              
              # Create sequence name
              seq_name <- paste0(accession, "_", gene_name, "_", feature_info, "_", j)
              id_extracted_sequences[[seq_name]] <- gene_subseq
              
              cat("Successfully extracted", gene_name, "(", length(gene_subseq), "bp)\n")
            } else {
              cat("WARNING: Invalid coordinates for", gene_name, "- skipping\n")
            }
          }
        }, error = function(e) {
          cat("ERROR: Failed to extract", gene_name, ":", e$message, "\n")
        })
      }
      
      # Save extracted sequences if any were found
      if (length(id_extracted_sequences) > 0) {
        cat("Saving", length(id_extracted_sequences), "extracted sequences for accession:", accession, "\n")
        
        gene_stringset <- Biostrings::DNAStringSet(id_extracted_sequences)
        
        # Create clean FASTA headers: accession|gene_name|[sequence_type]
        seq_type <- sequence_info[[id]]$type
        clean_names <- character(length(id_extracted_sequences))
        
        for (k in seq_along(id_extracted_sequences)) {
          seq_name <- names(id_extracted_sequences)[k]
          # Extract gene name from sequence name (format: accession_genename_featuretype_number)
          parts <- strsplit(seq_name, "_")[[1]]
          gene_name <- parts[2]  # Second part is the gene name
          
          # Create clean header format
          clean_names[k] <- paste0(accession, "|", gene_name, "|[", seq_type, "]")
        }
        names(gene_stringset) <- clean_names
        
        # Create output filename
        output_filename <- paste0(taxon_name_clean, "_taxid", taxid, "_", accession, "_extracted_genes.fasta")
        output_file <- file.path(extracted_folder, output_filename)
        
        # Write sequences to file
        Biostrings::writeXStringSet(gene_stringset, output_file)
        files_created <- files_created + 1
        
        cat("Saved sequences to:", output_file, "\n")
      } else {
        cat("No sequences could be extracted for accession:", accession, "\n")
      }
    }
  }
  
  cat("Created", files_created, "FASTA files with extracted sequences\n")
  return(files_created)
}

# =============================================================================
# MAIN PROCESSING FUNCTIONS
# =============================================================================

#' Process a single taxonomic ID
#'
#' @param taxid Character string containing taxonomic ID
#' @param gene_synonyms Character vector of gene names to search for
#' @param feature_type Character string specifying feature type filter
#' @param retmax Integer maximum number of sequences to retrieve
#' @param output_base_dir Base directory for output files
#' @return List containing processing results for this taxid
process_single_taxid <- function(taxid, gene_synonyms, feature_type = "all", 
                                 retmax = 5, output_base_dir) {
  
  cat("\n", paste(rep("=", 50), collapse = ""), "\n")
  cat("PROCESSING TAXID:", taxid, "\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  
  tryCatch({
    # Get taxonomic information
    taxon_info <- get_taxon_info(taxid)
    taxon_name <- taxon_info$taxon_name
    taxon_name_clean <- taxon_info$taxon_name_clean
    
    cat("Taxon name:", taxon_name, "\n")
    
    # Search NCBI for sequences
    search_results <- search_ncbi_sequences(taxid, gene_synonyms, retmax)
    ids <- search_results$ids
    full_query <- search_results$query
    
    if (length(ids) == 0) {
      cat("No sequences found for taxid", taxid, "\n")
      return(list(
        taxid = taxid,
        taxon_name = taxon_name,
        query = full_query,
        ids = character(0),
        message = "No sequences found"
      ))
    }
    
    cat("Found", length(ids), "sequences for taxid", taxid, "\n")
    
    # Create output directories
    dirs <- create_output_directories(output_base_dir)
    
    # Download and process sequences
    seq_data <- download_and_process_sequences(ids, dirs$raw_files_folder)
    
    # Create sequence type summary
    seq_summary_df <- do.call(rbind, lapply(seq_data$sequence_info, function(x) {
      data.frame(
        Taxid = taxid,
        Taxon_Name = taxon_name,
        Sequence_ID = x$id,
        Accession = x$accession,
        Sequence_Type = x$type,
        Title = x$title,
        stringsAsFactors = FALSE
      )
    }))
    
    # Print sequence type summary
    if (nrow(seq_summary_df) > 0) {
      cat("\nSequence Type Summary for Taxid", taxid, ":\n")
      type_counts <- table(seq_summary_df$Sequence_Type)
      for (type in names(type_counts)) {
        cat("  ", type, ":", type_counts[type], "\n")
      }
    }
    
    # Filter gene coordinates
    gene_coords <- filter_gene_coordinates(seq_data$gb_str_dfs, gene_synonyms, feature_type)
    
    # Extract and save gene sequences
    files_created <- extract_and_save_gene_sequences(
      gene_coords, seq_data$gb_seq_dfs, seq_data$sequence_info,
      seq_data$id_to_accession, dirs$extracted_folder,
      taxon_name_clean, taxid
    )
    
    # Create analysis summary file
    create_analysis_summary(taxid, taxon_name, taxon_name_clean, full_query, 
                            feature_type, seq_data, gene_coords, files_created,
                            dirs$output_folder, type_counts)
    
    cat("Successfully completed processing taxid", taxid, "\n")
    
    # Return comprehensive results
    return(list(
      taxid = taxid,
      taxon_name = taxon_name,
      query = full_query,
      ids = ids,
      id_to_accession = seq_data$id_to_accession,
      gb_str_dfs = seq_data$gb_str_dfs,
      fasta_sequences = seq_data$gb_seq_dfs,
      sequence_info = seq_data$sequence_info,
      sequence_summary = seq_summary_df,
      gene_coordinates = gene_coords,
      output_folder = dirs$output_folder,
      files_created = files_created,
      type_counts = if(exists("type_counts")) type_counts else NULL
    ))
    
  }, error = function(e) {
    cat("ERROR: Failed to process taxid", taxid, ":", e$message, "\n")
    return(list(
      taxid = taxid,
      error = e$message
    ))
  })
}

#' Create analysis summary file for a processed taxid
#' 
#' @param taxid Taxonomic ID
#' @param taxon_name Scientific name of taxon
#' @param taxon_name_clean Cleaned name for file naming
#' @param full_query NCBI search query used
#' @param feature_type Feature type filter applied
#' @param seq_data Processed sequence data
#' @param gene_coords Filtered gene coordinates
#' @param files_created Number of output files created
#' @param output_folder Output directory path
#' @param type_counts Table of sequence type counts
create_analysis_summary <- function(taxid, taxon_name, taxon_name_clean, full_query,
                                    feature_type, seq_data, gene_coords, files_created,
                                    output_folder, type_counts) {
  
  summary_file <- file.path(output_folder, paste0(taxon_name_clean, "_taxid", taxid, "_analysis_summary.txt"))
  
  cat("Creating analysis summary file:", summary_file, "\n")
  
  # Write summary information
  cat("Analysis Summary for Taxid:", taxid, "\n", file = summary_file)
  cat("Taxon name:", taxon_name, "\n", file = summary_file, append = TRUE)
  cat("Query:", full_query, "\n", file = summary_file, append = TRUE)
  cat("Feature type filter:", feature_type, "\n", file = summary_file, append = TRUE)
  cat("Total sequences processed:", length(seq_data$sequence_info), "\n", file = summary_file, append = TRUE)
  cat("GenBank files created:", length(seq_data$gb_str_dfs), "\n", file = summary_file, append = TRUE)
  cat("FASTA files created:", length(seq_data$gb_seq_dfs), "\n", file = summary_file, append = TRUE)
  cat("Gene matches found:", length(gene_coords), "\n", file = summary_file, append = TRUE)
  cat("Separate FASTA files created:", files_created, "\n", file = summary_file, append = TRUE)
  
  # Add sequence type breakdown
  if (!is.null(type_counts) && length(type_counts) > 0) {
    cat("\nSequence Type Breakdown:\n", file = summary_file, append = TRUE)
    for (type in names(type_counts)) {
      cat("  ", type, ":", type_counts[type], "\n", file = summary_file, append = TRUE)
    }
  }
  
  # Add accession mapping
  cat("\nNCBI ID to Accession Mapping:\n", file = summary_file, append = TRUE)
  for (ncbi_id in names(seq_data$id_to_accession)) {
    cat("  ", ncbi_id, "->", seq_data$id_to_accession[[ncbi_id]], "\n", file = summary_file, append = TRUE)
  }
  
  cat("Analysis summary saved to:", summary_file, "\n")
}

#' Print overall summary of processing results
#' 
#' @param all_results List of results from all processed taxids
#' @param overall_sequence_summary Data frame containing all sequence information
#' @param output_base_dir Base output directory
print_overall_summary <- function(all_results, overall_sequence_summary, output_base_dir) {
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("OVERALL PROCESSING SUMMARY\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  # Calculate success/failure statistics
  total_taxids <- length(all_results)
  successful <- sum(sapply(all_results, function(x) is.null(x$error)))
  failed <- sum(sapply(all_results, function(x) !is.null(x$error)))
  
  cat("Total taxids processed:", total_taxids, "\n")
  cat("Successful:", successful, "\n")
  cat("Failed:", failed, "\n")
  
  # Print details of failed taxids
  if (failed > 0) {
    cat("\nFailed taxids:\n")
    for (taxid in names(all_results)) {
      result <- all_results[[taxid]]
      if (!is.null(result$error)) {
        cat("  ", taxid, ":", result$error, "\n")
      }
    }
  }
  
  # Print overall sequence type summary
  if (nrow(overall_sequence_summary) > 0) {
    cat("\nOVERALL SEQUENCE TYPE SUMMARY:\n")
    overall_type_counts <- table(overall_sequence_summary$Sequence_Type)
    for (type in names(overall_type_counts)) {
      cat("  ", type, ":", overall_type_counts[type], "\n")
    }
    
    # Save overall summary table
    overall_summary_file <- file.path(output_base_dir, "overall_sequence_summary.csv")
    write.csv(overall_sequence_summary, overall_summary_file, row.names = FALSE)
    cat("\nOverall sequence summary saved to:", overall_summary_file, "\n")
    
    if (nrow(overall_sequence_summary) > 10) {
      cat("... and", nrow(overall_sequence_summary) - 10, "more rows\n")
    }
  } else {
    cat("\nNo sequences were successfully processed.\n")
  }
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("WORKFLOW COMPLETED\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
}

#' Main function to extract gene sequences from multiple taxonomic IDs
#' 
#' @param taxid_file Path to file containing taxonomic IDs (one per line)
#' @param gene_synonyms Character vector of gene names/synonyms to search for
#' @param feature_type Character string specifying feature type filter ("all" for no filtering)
#' @param retmax Integer maximum number of sequences to retrieve per taxid
#' @param output_base_dir Base directory for all output files
#' @return Data frame containing overall summary of all processed sequences
#' @export
extract_gene_sequences <- function(taxid_file, gene_synonyms, feature_type = "all", 
                                   retmax = 5, output_base_dir) {
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("STARTING GENE SEQUENCE EXTRACTION WORKFLOW\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  cat("Parameters:\n")
  cat("  Taxid file:", taxid_file, "\n")
  cat("  Gene synonyms:", paste(gene_synonyms, collapse = ", "), "\n")
  cat("  Feature type filter:", feature_type, "\n")
  cat("  Max sequences per taxid:", retmax, "\n")
  cat("  Output directory:", output_base_dir, "\n\n")
  
  # Read and validate taxonomic IDs
  taxids <- read_taxids(taxid_file)
  
  # Initialize results storage
  all_results <- list()
  overall_sequence_summary <- data.frame()
  
  # Process each taxonomic ID
  for (i in seq_along(taxids)) {
    taxid <- taxids[i]
    cat("Processing taxid", taxid, "(", i, "of", length(taxids), ")\n")
    
    # Process single taxid
    result <- process_single_taxid(taxid, gene_synonyms, feature_type, retmax, output_base_dir)
    all_results[[taxid]] <- result
    
    # Add to overall summary if successful
    if (!is.null(result$sequence_summary) && nrow(result$sequence_summary) > 0) {
      overall_sequence_summary <- rbind(overall_sequence_summary, result$sequence_summary)
    }
  }
  
  # Print overall processing summary
  print_overall_summary(all_results, overall_sequence_summary, output_base_dir)
  
  # Return the overall sequence summary
  return(overall_sequence_summary)
}

# =============================================================================
# test
# =============================================================================
coi_synonyms <- c("COI[Gene]",
                  "COX1[Gene]",
                  "cytochrome c oxidase subunit 1[Title]",
                  "cytochrome oxidase[All Fields]")

extract_gene_sequences("taxids/proseriate_taxid.txt", coi_synonyms, "CDS", 10, "test/new_proseriate_test")
extract_gene_sequences("taxids/metazoans_5_taxids.txt", coi_synonyms, "CDS", 3, "test/new_metazoan_test")
