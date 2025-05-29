# Load required libraries
library(rentrez)  # For interfacing with NCBI Entrez databases

# Improved helper function to extract gene_target features from GenBank records
extract_gene_target_from_genbank <- function(genbank_text, genome_id, gene_synonyms) {
  
  # Split GenBank record into lines
  lines <- strsplit(genbank_text, "\n")[[1]]
  
  # Initialize variables
  gene_target_sequences <- character()
  in_features <- FALSE
  current_feature <- NULL
  sequence_lines <- character()
  in_origin <- FALSE
  gene_target_features <- list()
  
  # Parse GenBank record
  for (line in lines) {
    
    # Check if we're in the FEATURES section
    if (grepl("^FEATURES", line)) {
      in_features <- TRUE
      next
    }
    
    # Check if we're in the ORIGIN section (sequence data)
    if (grepl("^ORIGIN", line)) {
      in_origin <- TRUE
      in_features <- FALSE
      next
    }
    
    # End of record
    if (grepl("^//", line)) {
      break
    }
    
    # Collect sequence data
    if (in_origin) {
      # Remove line numbers and spaces, keep only nucleotides
      clean_line <- gsub("[^acgtACGTnN]", "", line)
      if (nchar(clean_line) > 0) {
        sequence_lines <- c(sequence_lines, clean_line)
      }
    }
    
    # Process features
    if (in_features && !in_origin) {
      
      # New feature line (starts with 5 spaces and feature type)
      if (grepl("^     [A-Za-z]", line)) {
        
        # Process previous feature if it was gene_target
        if (!is.null(current_feature) && is_gene_target_feature(current_feature$lines, gene_synonyms)) {
          location <- extract_feature_location(current_feature$lines)
          if (!is.null(location)) {
            current_feature$location <- location
            gene_target_features <- append(gene_target_features, list(current_feature))
          }
        }
        
        # Start new feature
        feature_type <- trimws(strsplit(line, "\\s+")[[1]][1])
        current_feature <- list(type = feature_type, lines = c(line))
        
      } else if (grepl("^                     ", line) && !is.null(current_feature)) {
        # Continuation of current feature (21 spaces)
        current_feature$lines <- c(current_feature$lines, line)
      }
    }
  }
  
  # Process the last feature
  if (!is.null(current_feature) && is_gene_target_feature(current_feature$lines, gene_synonyms)) {
    location <- extract_feature_location(current_feature$lines)
    if (!is.null(location)) {
      current_feature$location <- location
      gene_target_features <- append(gene_target_features, list(current_feature))
    }
  }
  
  # Reconstruct full genome sequence
  full_sequence <- paste(sequence_lines, collapse = "")
  
  # Extract gene_target sequences based on coordinates
  extracted_sequences <- character()
  
  if (length(gene_target_features) > 0) {
    for (i in seq_along(gene_target_features)) {
      feature <- gene_target_features[[i]]
      
      if (!is.null(feature$location)) {
        tryCatch({
          # Extract sequence based on coordinates
          start_pos <- feature$location$start
          end_pos <- feature$location$end
          
          if (start_pos <= nchar(full_sequence) && end_pos <= nchar(full_sequence) && start_pos > 0) {
            gene_target_seq <- substr(full_sequence, start_pos, end_pos)
            
            # Handle complement sequences if needed
            if (feature$location$complement) {
              gene_target_seq <- reverse_complement(gene_target_seq)
            }
            
            # Create FASTA header
            header <- paste0(">gene_target_from_", genome_id, "_", start_pos, "..", end_pos, 
                             ifelse(feature$location$complement, "_complement", ""))
            
            extracted_sequences <- c(extracted_sequences, header, gene_target_seq)
          }
          
        }, error = function(e) {
          cat("      Error extracting gene_target sequence:", e$message, "\n")
        })
      }
    }
  }
  
  return(extracted_sequences)
}

# Helper function to check if a feature is gene_target-related
is_gene_target_feature <- function(feature_lines, gene_synonyms) {
  # Convert all lines to a single string for searching
  feature_text <- paste(feature_lines, collapse = " ")
  
  # Look for gene_target names in various forms
  for (pattern in gene_synonyms) {
    if (grepl(pattern, feature_text, ignore.case = TRUE)) {
      return(TRUE)
    }
  }
  
  return(FALSE)
}

# Improved helper function to extract feature location coordinates
extract_feature_location <- function(feature_lines) {
  
  # Look for location information in the feature lines
  location_text <- paste(feature_lines, collapse = " ")
  
  # Initialize complement flag
  is_complement <- FALSE
  
  # Check for complement
  if (grepl("complement", location_text, ignore.case = TRUE)) {
    is_complement <- TRUE
  }
  
  # Extract coordinates from various patterns
  # Pattern 1: Simple coordinates like "570..9008"
  simple_match <- regexpr("\\d+\\.\\.\\d+", location_text)
  
  if (simple_match > 0) {
    location_str <- regmatches(location_text, simple_match)
    coords <- strsplit(location_str, "\\.\\.")[[1]]
    
    if (length(coords) == 2) {
      start_pos <- as.numeric(coords[1])
      end_pos <- as.numeric(coords[2])
      
      # Validate coordinates
      if (!is.na(start_pos) && !is.na(end_pos) && start_pos > 0 && end_pos > start_pos) {
        return(list(start = start_pos, end = end_pos, complement = is_complement))
      }
    }
  }
  
  return(NULL)
}

# Simple reverse complement function
reverse_complement <- function(seq) {
  # Convert to uppercase
  seq <- toupper(seq)
  
  # Create complement
  complement_map <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G", "N" = "N")
  seq_chars <- strsplit(seq, "")[[1]]
  complement_chars <- complement_map[seq_chars]
  complement_chars[is.na(complement_chars)] <- "N"  # Handle any unknown characters
  
  # Reverse the sequence
  rev_complement <- paste(rev(complement_chars), collapse = "")
  
  return(rev_complement)
}

create_gene_search_query <- function(taxid, gene_synonyms) {
  # Create OR conditions for each gene synonym
  gene_conditions <- character()
  
  for (synonym in gene_synonyms) {
    # Add different search patterns for each synonym
    gene_conditions <- c(gene_conditions, 
                         paste0(synonym, "[Gene]"),
                         paste0("\"", synonym, "\"[All Fields]"))
  }
  
  # Combine all conditions with OR
  gene_query_part <- paste(gene_conditions, collapse = " OR ")
  
  # Complete query
  full_query <- paste0("txid", taxid, "[Organism] AND (", gene_query_part, ")")
  
  return(full_query)
}


get_gene_target_by_order <- function(taxid_file, max_per_order, output_dir, gene_synonyms) {
  
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
  cat("Processing", length(taxids), "taxonomic orders\n")
  cat("Getting up to", max_per_order, "sequences per order\n")
  cat("Searching for genes:", paste(gene_synonyms, collapse = ", "), "\n")
  
  # Initialize results tracking dataframe
  results_summary <- data.frame(
    taxid = character(),
    order_name = character(),
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
      
      # Search for gene sequences
      search_result <- rentrez::entrez_search(
        db = "nuccore",
        term = search_query,
        retmax = max_per_order * 3  # Get more IDs to account for filtering
      )
      
      sequences_found <- search_result$count
      sequences_to_process <- length(search_result$ids)
      
      cat("  Found:", sequences_found, "total sequences\n")
      cat("  Will process:", sequences_to_process, "sequences\n")
      
      # Get taxonomic name
      order_name <- "Unknown"
      tax_summary <- rentrez::entrez_summary(db = "taxonomy", id = taxid)
      if (!is.null(tax_summary$scientificname)) {
        order_name <- tax_summary$scientificname
      }
      cat("  Order name:", order_name, "\n")
      
      # Initialize variables for this taxon
      final_sequences <- character()
      sequences_retrieved <- 0
      gene_target_extracted_count <- 0
      
      # Process sequences up to max_per_order
      if (sequences_to_process > 0) {
        for (seq_id in search_result$ids) {
          if (sequences_retrieved >= max_per_order) break
          
          tryCatch({
            # Get sequence summary to check if it's a complete genome
            seq_summary <- rentrez::entrez_summary(db = "nuccore", id = seq_id)
            
            # Check if this is a complete genome
            if (!is.null(seq_summary$title) && grepl("complete genome", seq_summary$title, ignore.case = TRUE)) {
              
              # This is a complete genome - extract gene features
              cat("    Processing complete genome:", seq_id, "\n")
              
              # Fetch GenBank record
              genbank_record <- rentrez::entrez_fetch(
                db = "nuccore",
                id = seq_id,
                rettype = "gb"
              )
              
              # Save GenBank record for reference
              genbank_file <- file.path(genbank_dir, paste0("genome_", seq_id, ".gb"))
              writeLines(genbank_record, genbank_file)
              
              # Extract gene sequences from this genome using gene synonyms
              gene_sequences <- extract_gene_target_from_genbank(genbank_record, seq_id, gene_synonyms)
              
              if (length(gene_sequences) > 0) {
                final_sequences <- c(final_sequences, gene_sequences)
                extracted_count <- length(gene_sequences) / 2  # Divide by 2 (header + sequence)
                gene_target_extracted_count <- gene_target_extracted_count + extracted_count
                sequences_retrieved <- sequences_retrieved + extracted_count
                cat("      Extracted", extracted_count, "gene features\n")
              }
              
            } else {
              
              # This is a regular gene sequence - fetch directly
              gene_sequence <- rentrez::entrez_fetch(
                db = "nuccore",
                id = seq_id,
                rettype = "fasta"
              )
              
              final_sequences <- c(final_sequences, gene_sequence)
              sequences_retrieved <- sequences_retrieved + 1
              cat("    Retrieved gene sequence:", seq_id, "\n")
            }
            
            # Small delay between requests
            Sys.sleep(0.3)
            
          }, error = function(e) {
            cat("    Error processing sequence", seq_id, ":", e$message, "\n")
          })
        }
      }
      
      # Save sequences if any were obtained
      if (length(final_sequences) > 0) {
        safe_order_name <- gsub("[^A-Za-z0-9]", "_", order_name)
        output_file <- file.path(output_dir, paste0("taxid_", taxid, "_", safe_order_name, ".fasta"))
        
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
        order_name = order_name,
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
        order_name = "ERROR",
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
  cat("Total orders processed:", nrow(results_summary), "\n")
  cat("Orders with sequences:", sum(results_summary$sequences_retrieved > 0), "\n")
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

# EXAMPLE USAGE:

# Define gene synonyms for COI/COX1
coi_synonyms <- c(
  "COI", "COX1", "COXI",
  "cytochrome c oxidase subunit 1",
  "cytochrome c oxidase subunit I",
  "cytochrome oxidase subunit 1",
  "cytochrome oxidase subunit I"
)

# Or define synonyms for a different gene (example: 16S rRNA)
# rRNA_16S_synonyms <- c(
#   "16S", "16S rRNA", "16S ribosomal RNA",
#   "small subunit ribosomal RNA",
#   "SSU rRNA"
# )

# Run the analysis with COI synonyms
results <- get_gene_target_by_order(
  taxid_file = "all_animal_order_taxids.txt",
  max_per_order = 2,
  output_dir = "all_animal_orders",
  gene_synonyms = coi_synonyms  # Pass the gene synonyms vector here
)

# Combine all sequences into one file
combine_sequences(
  output_dir = "all_animal_orders",
  combined_file = "all_animal_orders/gene_sequences.fasta"
)