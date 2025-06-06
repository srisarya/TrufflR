# TrufflR: Extract specific gene sequences from NCBI's database for specific taxa

# USER INPUTS:
# - taxid_file: Text file containing NCBI taxonomic IDs (one per line)
# - gene_synonyms: Vector of search terms with NCBI field tags (e.g., "COI[Gene]")
# - max_per_taxid: Maximum sequences to retrieve per taxonomic group
# - output_dir: Directory path for saving results
# - feature_types: Vector of GenBank feature types to search (e.g., c("CDS", "rRNA"))

# OUTPUTS:
# - Individual FASTA files per taxonomic group
# - Combined master FASTA file with all sequences
# - CSV summary with retrieval statistics
# - GenBank records folder with complete genome files
# - Error log for troubleshooting failed retrievals


# Load required libraries
library(rentrez)     # For interfacing with NCBI Entrez databases
library(geneviewer)  # For parsing GenBank files
library(Biostrings)  # For sequence manipulation
library(seqinr)      # For sequence manipulation and reverse complement

# Functions
create_gene_search_query <- function(taxid, gene_synonyms) {
  # Use gene synonyms exactly as provided by user (no additional field tags)
  # Combine all conditions with OR
  gene_query_part <- paste(gene_synonyms, collapse = " OR ")
  
  # Complete query - only add the taxid part
  full_query <- paste0("txid", taxid, "[Organism] AND (", gene_query_part, ")")
  
  return(full_query)
}

# Function to validate if a sequence record matches target gene synonyms
validate_gene_match <- function(seq_summary, gene_synonyms) {
  # Clean gene synonyms by removing NCBI field tags
  clean_synonyms <- tolower(gsub("\\[.*?\\]", "", gene_synonyms))
  
  # Get title and other relevant fields from summary
  title <- tolower(seq_summary$title %||% "")
  
  # Check if any synonym appears in the title
  for (synonym in clean_synonyms) {
    if (grepl(synonym, title, ignore.case = TRUE)) {
      return(TRUE)
    }
  }
  
  return(FALSE)
}

# Function to determine sequence type based on title
determine_sequence_type <- function(title) {
  title_lower <- tolower(title)
  
  if (grepl("partial", title_lower)) {
    return("PARTIAL_CDS")
  } else if (grepl("genome|chromosome|scaffold", title_lower)) {
    return("GENOMIC")
  } else {
    return("COMPLETE_CDS")
  }
}

# Updated function using geneviewer's gbk_features_to_df for feature extraction
extract_gene_target_with_geneviewer <- function(genbank_file, genome_id, gene_synonyms, feature_types = c("CDS"), return_aa = FALSE) {
  tryCatch({
    # Read GenBank file using geneviewer - handle both file path and temp file
    if (file.exists(genbank_file)) {
      gb_data <- geneviewer::read_gbk(genbank_file)
    } else {
      stop("GenBank file not found")
    }
    
    # Initialize vector to store all extracted sequences
    all_extracted_sequences <- character()
    
    # Process each feature type
    for (feature_type in feature_types) {
      # Extract features of specified type using gbk_features_to_df
      features_df <- geneviewer::gbk_features_to_df(
        genbank_file,
        feature = feature_type,
        keys = NULL,
        process_region = TRUE
      )
      
      if (nrow(features_df) == 0) {
        cat("    No", feature_type, "features found in genome\n")
        next
      }
      
      cat("    Feature columns available:", paste(colnames(features_df), collapse = ", "), "\n")
      
      # Prepare gene synonyms for matching (remove NCBI field tags and convert to lowercase)
      clean_synonyms <- tolower(gsub("\\[.*?\\]", "", gene_synonyms))
      cat("    Looking for genes:", paste(clean_synonyms, collapse = ", "), "\n")
      
      # Look for matching gene names in available columns
      gene_matches <- c()
      
      # Check 'gene' column if it exists
      if ("gene" %in% colnames(features_df)) {
        gene_col_matches <- which(sapply(features_df$gene, function(x) {
          if (is.na(x)) return(FALSE)
          any(sapply(clean_synonyms, function(syn) grepl(syn, tolower(x), ignore.case = TRUE)))
        }))
        gene_matches <- c(gene_matches, gene_col_matches)
        if (length(gene_col_matches) > 0) {
          cat("    Found matches in 'gene' column:", gene_col_matches, "\n")
        }
      }
      
      # Check 'product' column if it exists
      if ("product" %in% colnames(features_df)) {
        product_col_matches <- which(sapply(features_df$product, function(x) {
          if (is.na(x)) return(FALSE)
          any(sapply(clean_synonyms, function(syn) grepl(syn, tolower(x), ignore.case = TRUE)))
        }))
        gene_matches <- c(gene_matches, product_col_matches)
        if (length(product_col_matches) > 0) {
          cat("    Debug: Found matches in 'product' column:", product_col_matches, "\n")
        }
      }
      
      # Remove duplicates
      all_matches <- unique(gene_matches)
      
      if (length(all_matches) == 0) {
        cat("    No matching genes found in", feature_type, "features\n")
        next
      }
      
      cat("    Found", length(all_matches), "matching", feature_type, "features\n")
      
      # Get the complete sequence for extraction
      complete_seq <- NULL
      if (!is.null(gb_data$sequence) && nchar(gb_data$sequence) > 0) {
        complete_seq <- gb_data$sequence
      } else {
        warning("No sequence data found in GenBank record")
        next
      }
      
      # Extract sequences for all matches
      for (match_idx in all_matches) {
        feature <- features_df[match_idx, ]
        cat("    Debug: Processing feature at row", match_idx, "\n")
        
        # For CDS features, check if we want amino acid sequence
        if (feature_type == "CDS" && return_aa && "translation" %in% colnames(features_df)) {
          # Use translated amino acid sequence if available
          gene_seq <- feature$translation
          if (is.na(gene_seq) || gene_seq == "") {
            cat("    Warning: No translation found for CDS feature, falling back to DNA\n")
            gene_seq <- NULL
          } else {
            cat("    Debug: Using amino acid translation\n")
          }
        } else {
          gene_seq <- NULL
        }
        
        # If no amino acid sequence or not requested, extract DNA sequence
        if (is.null(gene_seq)) {
          # Get sequence coordinates from the dataframe
          start_pos <- feature$start
          end_pos <- feature$end
          strand <- feature$strand
          
          cat("    Debug: Extracting DNA from positions", start_pos, "to", end_pos, "strand:", strand, "\n")
          
          # Validate coordinates
          if (is.na(start_pos) || is.na(end_pos) || start_pos < 1 || end_pos > nchar(complete_seq)) {
            cat("    Warning: Invalid coordinates for feature - start:", start_pos, "end:", end_pos, "genome length:", nchar(complete_seq), "\n")
            next
          }
          
          # Extract subsequence
          gene_seq <- substr(complete_seq, start_pos, end_pos)
          
          # Apply reverse complement if on negative strand using seqinr
          if (!is.na(strand) && (strand == "-" || strand == -1)) {
            gene_seq <- seqinr::c2s(seqinr::comp(seqinr::s2c(gene_seq), forceToLower = FALSE))
            gene_seq <- seqinr::c2s(rev(seqinr::s2c(gene_seq)))
            cat("    Debug: Applied reverse complement using seqinr\n")
          }
        }
        
        # Create FASTA header with feature information
        gene_name <- "unknown"
        if ("gene" %in% colnames(features_df) && !is.na(feature$gene) && feature$gene != "") {
          gene_name <- feature$gene
        } else if ("product" %in% colnames(features_df) && !is.na(feature$product) && feature$product != "") {
          gene_name <- feature$product
        }
        
        seq_type <- if (feature_type == "CDS" && return_aa && "translation" %in% colnames(features_df) && !is.na(feature$translation)) "AA" else "DNA"
        header <- paste0(">", genome_id, "_taxid_GENOME_EXTRACTED_", gene_name, "_", feature_type, "_", seq_type, "_", 
                         feature$start, "..", feature$end, 
                         ifelse(!is.na(feature$strand) && (feature$strand == "-" || feature$strand == -1), "_complement", ""))
        
        all_extracted_sequences <- c(all_extracted_sequences, header, gene_seq)
        cat("    Debug: Added sequence of length", nchar(gene_seq), "\n")
      }
    }
    
    if (length(all_extracted_sequences) == 0) {
      stop("No matching genes found in any specified feature types")
    }
    
    return(all_extracted_sequences)
    
  }, error = function(e) {
    stop("geneviewer parsing failed: ", e$message)
  })
}

# Main function for gene extraction by taxonomic ID
get_gene_target_by_taxid <- function(taxid_file, max_per_taxid, output_dir, gene_synonyms, feature_types = c("CDS")) {
  
  # Validate gene_synonyms parameter
  if (missing(gene_synonyms) || length(gene_synonyms) == 0) {
    stop("gene_synonyms parameter is required and must contain at least one gene name")
  }
  
  # Read taxonomy IDs from input file
  if (!file.exists(taxid_file)) {
    stop("File ", taxid_file, " not found!")
  }
  
  # Read all lines from the taxonomy ID file
  taxids <- readLines(taxid_file)
  taxids <- taxids[taxids != ""]  # Remove empty lines
  
  # Print initial information about the analysis
  cat("Processing", length(taxids), "taxonomic IDs\n")
  cat("Getting up to", max_per_taxid, "sequences per taxid\n")
  cat("Searching for genes:", paste(gene_synonyms, collapse = ", "), "\n")
  cat("Looking for feature types:", paste(feature_types, collapse = ", "), "\n")
  
  # Initialize results tracking dataframe
  results_summary <- data.frame(
    taxid = character(),
    taxon_name = character(),
    sequences_found = integer(),
    sequences_retrieved = integer(),
    gene_target_extracted_from_genomes = integer(),
    stringsAsFactors = FALSE
  )
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Create subdirectory for GenBank records
  genbank_dir <- file.path(output_dir, "genbank_records")
  if (!dir.exists(genbank_dir)) {
    dir.create(genbank_dir)
  }
  
  # Process each taxonomy ID sequentially
  for (i in seq_along(taxids)) {
    taxid <- trimws(taxids[i])
    cat("Processing taxid", taxid, "(", i, "of", length(taxids), ")\n")
    
    tryCatch({
      
      # Construct search query using gene synonyms
      search_query <- create_gene_search_query(taxid, gene_synonyms)
      
      cat("  Query:", search_query, "\n")
      
      # Search for gene sequences with up to 3 retries
      search_result <- NULL
      attempt <- 1
      taxon_name <- "Unknown"
      
      while (is.null(search_result) && attempt <= 3) {
        tryCatch({
          search_result <- rentrez::entrez_search(
            db = "nuccore",
            term = search_query,
            retmax = max_per_taxid * 3  # Get more IDs to account for filtering
          )
          # Get taxonomic name
          tax_summary <- rentrez::entrez_summary(db = "taxonomy", id = taxid)
          if (!is.null(tax_summary$scientificname)) {
            taxon_name <- tax_summary$scientificname
          }
          cat("  Taxon name:", taxon_name, "\n")
        }, error = function(e) {
          cat("  Search attempt", attempt, "failed:", e$message, "\n")
          Sys.sleep(2 * attempt)  # Exponential backoff
        })
        attempt <- attempt + 1
      }
      
      if (is.null(search_result)) {
        write(paste("Failed to retrieve search results for taxid:", taxid, "taxon:", taxon_name), 
              file = file.path(output_dir, "error.log"), append = TRUE)
        stop("Failed to retrieve search results after 3 attempts")
      }
      
      sequences_found <- search_result$count
      sequences_to_process <- length(search_result$ids)
      
      cat("  Found:", sequences_found, "total sequences\n")
      cat("  Will process:", sequences_to_process, "sequences\n")
      
      # Initialize variables for this taxon
      final_sequences <- character()
      sequences_retrieved <- 0
      gene_target_extracted_count <- 0
      
      # Process sequences up to max_per_taxid
      if (sequences_to_process > 0) {
        for (seq_id in search_result$ids) {
          if (sequences_retrieved >= max_per_taxid) break
          
          attempt_seq <- 1
          success_seq <- FALSE
          
          while (attempt_seq <= 3 && !success_seq) {
            tryCatch({
              # Get sequence summary to check if it's a complete genome
              seq_summary <- rentrez::entrez_summary(db = "nuccore", id = seq_id)
              
              # Determine sequence type
              seq_type <- determine_sequence_type(seq_summary$title)
              
              # Check if this is a complete genome
              if (seq_type == "COMPLETE_GENOME") {
                
                # This is a complete genome - extract gene features using geneviewer
                cat("    Processing complete genome:", seq_id, "\n")
                
                # Fetch GenBank record with features
                genbank_record <- rentrez::entrez_fetch(
                  db = "nuccore",
                  id = seq_id,
                  rettype = "gbwithparts",
                  retmode = "text"
                )
                
                # Save GenBank record as temporary file for geneviewer
                temp_gbk <- tempfile(fileext = ".gbk")
                writeLines(genbank_record, temp_gbk)
                
                # Extract gene sequences using geneviewer
                gene_sequences <- extract_gene_target_with_geneviewer(temp_gbk, seq_id, gene_synonyms, feature_types)
                
                if (length(gene_sequences) > 0) {
                  final_sequences <- c(final_sequences, gene_sequences)
                  extracted_count <- length(gene_sequences) / 2  # Divide by 2 (header + sequence)
                  gene_target_extracted_count <- gene_target_extracted_count + extracted_count
                  sequences_retrieved <- sequences_retrieved + extracted_count
                  cat("      Extracted", extracted_count, "gene features\n")
                }
                
                # Clean up temp file
                unlink(temp_gbk)
                
              } else {
                
                # This is a regular gene sequence - validate it matches our target genes
                cat("    Validating", seq_type, "sequence:", seq_id, "\n")
                
                # Check if the sequence summary actually matches our gene synonyms
                if (validate_gene_match(seq_summary, gene_synonyms)) {
                  cat("    ✓ Sequence matches target genes, fetching...\n")
                  
                  # Fetch the sequence
                  gene_sequence <- rentrez::entrez_fetch(
                    db = "nuccore",
                    id = seq_id,
                    rettype = "fasta"
                  )
                  
                  # Parse fasta header and sequence
                  lines <- strsplit(gene_sequence, "\n")[[1]]
                  header <- lines[1]
                  sequence <- paste(lines[-1], collapse = "")
                  
                  # Create new header with gene info, taxid, and sequence type
                  new_header <- paste0(">", seq_id, "_taxid_", taxid, "_", seq_type, "_", gsub("^>", "", header))
                  
                  # Combine header and sequence
                  final_sequences <- c(final_sequences, new_header, sequence)
                  sequences_retrieved <- sequences_retrieved + 1
                  cat("    Retrieved validated", seq_type, "sequence:", seq_id, "\n")
                  
                } else {
                  cat("    ✗ Sequence does not match target genes, skipping\n")
                }
              }
              
              # Small delay between requests
              Sys.sleep(0.3)
              success_seq <- TRUE
              
            }, error = function(e) {
              cat("    Attempt", attempt_seq, "failed for sequence", seq_id, ":", e$message, "\n")
              Sys.sleep(0.5 * attempt_seq)
              attempt_seq <<- attempt_seq + 1
              if (attempt_seq > 3) {
                cat("    Failed to process sequence", seq_id, "after 3 attempts\n")
                write(paste("Failed to process sequence:", seq_id, "after 3 attempts in taxon", taxon_name), 
                      file = file.path(output_dir, "error.log"), append = TRUE)
              }
            })
          }
        }
      }
      
      # Save sequences if any were obtained
      if (length(final_sequences) > 0) {
        safe_taxon_name <- gsub("[^A-Za-z0-9]", "_", taxon_name)
        output_file <- file.path(output_dir, paste0("taxid_", taxid, "_", safe_taxon_name, ".fasta"))
        
        writeLines(final_sequences, output_file)
        cat("  Retrieved", sequences_retrieved, "sequences total\n")
        cat("  Gene features extracted from genomes:", gene_target_extracted_count, "\n")
        cat("  Saved to:", output_file, "\n")
      } else {
        cat("  No sequences obtained\n")
      }
      
      # Add results to summary
      results_summary <- rbind(results_summary, data.frame(
        taxid = taxid,
        taxon_name = taxon_name,
        sequences_found = sequences_found,
        sequences_retrieved = sequences_retrieved,
        gene_target_extracted_from_genomes = gene_target_extracted_count,
        stringsAsFactors = FALSE
      ))
      
      cat("  ✓ Complete\n\n")
      Sys.sleep(1.0)
      
    }, error = function(e) {
      cat("  ✗ Error:", e$message, "\n\n")
      
      results_summary <<- rbind(results_summary, data.frame(
        taxid = taxid,
        taxon_name = "ERROR",
        sequences_found = 0,
        sequences_retrieved = 0,
        gene_target_extracted_from_genomes = 0,
        stringsAsFactors = FALSE
      ))
    })
  }
  
  # Save summary
  write.csv(results_summary, file.path(output_dir, "retrieval_summary.csv"), row.names = FALSE)
  
  # Print final summary
  cat("=== FINAL SUMMARY ===\n")
  cat("Total taxa processed:", nrow(results_summary), "\n")
  cat("Taxa with sequences:", sum(results_summary$sequences_retrieved > 0), "\n")
  cat("Total sequences retrieved:", sum(results_summary$sequences_retrieved), "\n")
  cat("Gene features extracted from genomes:", sum(results_summary$gene_target_extracted_from_genomes), "\n")
  
  cat("\nFiles saved in:", output_dir, "/\n")
  cat("GenBank records saved in:", genbank_dir, "/\n")
  cat("Summary saved as: retrieval_summary.csv\n")
  
  return(results_summary)
}

# Function to combine all individual FASTA files into one master file
combine_sequences <- function(output_dir, combined_file) {
  
  # Find all FASTA files in the output directory
  fasta_files <- list.files(output_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  # Check if any FASTA files were found
  if (length(fasta_files) == 0) {
    cat("No FASTA files found in", output_dir, "\n")
    return()
  }
  
  cat("Combining", length(fasta_files), "FASTA files into", combined_file, "\n")
  
  # Initialize vector to store all sequences
  all_sequences <- character()
  
  # Read each FASTA file and combine sequences
  for (file in fasta_files) {
    # Read all lines from current FASTA file
    sequences <- readLines(file)
    # Append to master sequence collection
    all_sequences <- c(all_sequences, sequences)
  }
  
  # Write combined sequences to master file
  writeLines(all_sequences, combined_file)
  
  # Count total number of sequences (FASTA headers start with ">")
  total_seqs <- length(grep("^>", all_sequences))
  cat("Combined file contains", total_seqs, "sequences\n")
  cat("Saved as:", combined_file, "\n")
}

# Example Usage
coi_synonyms <- c(
  "COI[Gene]",
  "COX1[Gene]", 
  "cytochrome c oxidase subunit 1[Title]",
  "cytochrome oxidase[All Fields]"
)

# # Example for rRNA searches
# rrna_16s_synonyms <- c(
#   "16S[Title]", 
#   "16S rRNA[Title]", 
#   "16S ribosomal RNA[rRNA]",
#   "small subunit ribosomal RNA[Title]",
#   "SSU rRNA[Title]"
# )

# Run the analysis with COI synonyms
results <- get_gene_target_by_taxid(
  taxid_file = "proseriate_taxid.txt",
  max_per_taxid = 50,
  output_dir = "test_proseriate_CO1",
  gene_synonyms = coi_synonyms,
  feature_types = c("CDS")  # Default for protein-coding genes
)

# Combine all sequences into one file
combine_sequences(
  output_dir = "test_proseriate_CO1",
  combined_file = "test_proseriate_CO1/gene_sequences.fasta"
)

# Example for rRNA extraction
results_rrna <- get_gene_target_by_taxid(
  taxid_file = "metazoans_5_taxids.txt",
  max_per_taxid = 10,
  output_dir = "test_metazoans_5_taxids_CO1",
  gene_synonyms = coi_synonyms,
  feature_types = c("CDS")  # Specify rRNA features
)

# Combine all sequences into one file
combine_sequences(
  output_dir = "test_metazoans_5_taxids_CO1",
  combined_file = "test_metazoans_5_taxids_CO1/gene_sequences.fasta"
)