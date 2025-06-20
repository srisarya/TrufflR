# TrufflR: Extract specific gene sequences from NCBI's database for specific taxa
# FIXED VERSION - Properly extracts gene sequences from genomes using coordinates

# Load required libraries
library(rentrez)     # For interfacing with NCBI Entrez databases
library(geneviewer)  # For parsing GenBank files
library(Biostrings)  # For sequence manipulation
library(seqinr)      # For sequence manipulation and reverse complement

# Functions
create_gene_search_query <- function(taxid, gene_synonyms) {
  gene_query_part <- paste(gene_synonyms, collapse = " OR ")
  full_query <- paste0("txid", taxid, "[Organism] AND (", gene_query_part, ")")
  return(full_query)
}

# Function to validate if a sequence record matches target gene synonyms
validate_gene_match <- function(seq_summary, gene_synonyms) {
  clean_synonyms <- tolower(gsub("\\[.*?\\]", "", gene_synonyms))
  title <- tolower(seq_summary$title %||% "")
  
  for (synonym in clean_synonyms) {
    if (grepl(synonym, title, ignore.case = TRUE)) {
      return(TRUE)
    }
  }
  return(FALSE)
}

# FIXED: More stringent genomic validation
validate_genomic_gene_match <- function(genbank_record_text, gene_synonyms) {
  clean_synonyms <- tolower(gsub("\\[.*?\\]", "", gene_synonyms))
  genbank_lower <- tolower(genbank_record_text)
  
  # Look for actual gene features (more specific patterns)
  for (synonym in clean_synonyms) {
    # Look for gene features like:
    #     gene            123..456
    #                     /gene="COI"
    # or  CDS             123..456
    #                     /gene="COI"
    gene_pattern <- paste0('/gene="[^"]*', synonym, '[^"]*"')
    product_pattern <- paste0('/product="[^"]*', synonym, '[^"]*"')
    
    if (grepl(gene_pattern, genbank_lower, ignore.case = TRUE) || 
        grepl(product_pattern, genbank_lower, ignore.case = TRUE)) {
      cat("    ✓ Found target gene '", synonym, "' in GenBank features\n", sep = "")
      return(TRUE)
    }
  }
  
  cat("    ✗ Target genes not found in GenBank features\n")
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

# COMPLETELY REWRITTEN: Proper gene extraction from genomes
extract_gene_target_direct_parsing <- function(gbk_filename, genome_id, gene_synonyms, feature_types) {
  tryCatch({
    
    cat("    Extracting gene sequences from GenBank file using direct parsing\n")
    
    # Load required libraries
    if (!requireNamespace("Biostrings", quietly = TRUE)) {
      stop("Biostrings package is required but not installed")
    }
    
    # Read the entire GenBank file
    genbank_content <- readLines(gbk_filename)
    cat("    Read GenBank file with", length(genbank_content), "lines\n")
    
    # EXTRACT ORIGIN SEQUENCE
    origin_start <- which(grepl("^ORIGIN", genbank_content))
    if (length(origin_start) == 0) {
      stop("No ORIGIN section found in GenBank file")
    }
    
    origin_lines <- genbank_content[(origin_start + 1):length(genbank_content)]
    end_marker <- which(grepl("^//", origin_lines))[1]
    if (!is.na(end_marker)) {
      origin_lines <- origin_lines[1:(end_marker - 1)]
    }
    
    # Extract sequence
    sequence_parts <- character()
    for (line in origin_lines) {
      clean_line <- gsub("^\\s*\\d+\\s*", "", line)
      clean_seq <- toupper(gsub("\\s", "", clean_line))
      if (nchar(clean_seq) > 0 && grepl("^[ATCGN]+$", clean_seq)) {
        sequence_parts <- c(sequence_parts, clean_seq)
      }
    }
    
    complete_seq <- paste(sequence_parts, collapse = "")
    genome_seq <- Biostrings::DNAString(complete_seq)
    cat("    Extracted complete sequence of length:", length(genome_seq), "\n")
    
    # PARSE FEATURES DIRECTLY FROM GENBANK
    cat("    Parsing features directly from GenBank content\n")
    
    # Find the FEATURES section
    features_start <- which(grepl("^FEATURES", genbank_content))
    if (length(features_start) == 0) {
      stop("No FEATURES section found in GenBank file")
    }
    
    # Get features section (from FEATURES to ORIGIN)
    features_end <- origin_start - 1
    features_lines <- genbank_content[features_start:features_end]
    
    cat("    Found FEATURES section with", length(features_lines), "lines\n")
    
    # Parse each feature type
    all_extracted_sequences <- character()
    
    for (feature_type in feature_types) {
      cat("    Processing feature type:", feature_type, "\n")
      
      # Find all lines that start with the feature type
      feature_pattern <- paste0("^\\s{5}", feature_type, "\\s+")
      feature_starts <- which(grepl(feature_pattern, features_lines))
      
      cat("    Found", length(feature_starts), "potential", feature_type, "features\n")
      
      if (length(feature_starts) == 0) {
        next
      }
      
      # Process each feature
      for (i in seq_along(feature_starts)) {
        start_line <- feature_starts[i]
        
        # Find the end of this feature (next feature or end of section)
        if (i < length(feature_starts)) {
          end_line <- feature_starts[i + 1] - 1
        } else {
          # Find next feature that starts at column 6 or end of features
          next_feature <- which(grepl("^\\s{5}\\S", features_lines[(start_line + 1):length(features_lines)]))
          if (length(next_feature) > 0) {
            end_line <- start_line + next_feature[1] - 1
          } else {
            end_line <- length(features_lines)
          }
        }
        
        # Extract this feature's lines
        feature_lines <- features_lines[start_line:end_line]
        
        cat("    DEBUG - Processing feature", i, "from lines", start_line, "to", end_line, "\n")
        cat("    DEBUG - First line:", feature_lines[1], "\n")
        
        # Parse coordinates from the first line
        coord_line <- feature_lines[1]
        
        # Extract coordinates - handle various formats
        # Examples: "CDS             123..456", "CDS             complement(123..456)"
        coord_match <- regexpr("\\d+\\.\\.\\d+", coord_line)
        if (coord_match == -1) {
          cat("    Could not parse coordinates from:", coord_line, "\n")
          next
        }
        
        coords <- regmatches(coord_line, coord_match)
        coord_parts <- strsplit(coords, "\\.\\.")[[1]]
        start_pos <- as.numeric(coord_parts[1])
        end_pos <- as.numeric(coord_parts[2])
        
        # Check if it's on complement strand
        is_complement <- grepl("complement", coord_line)
        
        cat("    DEBUG - Parsed coordinates:", start_pos, "to", end_pos, 
            ifelse(is_complement, "(complement)", ""), "\n")
        
        # Extract gene information from qualifier lines
        gene_info <- list(gene = "", product = "", locus_tag = "", note = "")
        
        for (line in feature_lines[-1]) {  # Skip the coordinate line
          line <- trimws(line)
          
          # Parse /gene= qualifier
          if (grepl("^/gene=", line)) {
            gene_match <- regexpr('/gene="([^"]*)"', line, perl = TRUE)
            if (gene_match != -1) {
              gene_info$gene <- regmatches(line, gene_match)
              gene_info$gene <- gsub('/gene="([^"]*)"', '\\1', gene_info$gene, perl = TRUE)
            }
          }
          
          # Parse /product= qualifier
          if (grepl("^/product=", line)) {
            product_match <- regexpr('/product="([^"]*)"', line, perl = TRUE)
            if (product_match != -1) {
              gene_info$product <- regmatches(line, product_match)
              gene_info$product <- gsub('/product="([^"]*)"', '\\1', gene_info$product, perl = TRUE)
            }
          }
          
          # Parse /locus_tag= qualifier
          if (grepl("^/locus_tag=", line)) {
            locus_match <- regexpr('/locus_tag="([^"]*)"', line, perl = TRUE)
            if (locus_match != -1) {
              gene_info$locus_tag <- regmatches(line, locus_match)
              gene_info$locus_tag <- gsub('/locus_tag="([^"]*)"', '\\1', gene_info$locus_tag, perl = TRUE)
            }
          }
          
          # Parse /note= qualifier
          if (grepl("^/note=", line)) {
            note_match <- regexpr('/note="([^"]*)"', line, perl = TRUE)
            if (note_match != -1) {
              gene_info$note <- regmatches(line, note_match)
              gene_info$note <- gsub('/note="([^"]*)"', '\\1', gene_info$note, perl = TRUE)
            }
          }
        }
        
        cat("    DEBUG - Gene info extracted:\n")
        cat("      gene:", gene_info$gene, "\n")
        cat("      product:", gene_info$product, "\n")
        cat("      locus_tag:", gene_info$locus_tag, "\n")
        cat("      note:", gene_info$note, "\n")
        
        # Check if this gene matches our search terms
        clean_synonyms <- tolower(gsub("\\[.*?\\]", "", gene_synonyms))
        match_found <- FALSE
        match_field <- ""
        match_value <- ""
        
        for (field in names(gene_info)) {
          if (gene_info[[field]] != "") {
            for (synonym in clean_synonyms) {
              if (grepl(synonym, tolower(gene_info[[field]]), ignore.case = TRUE)) {
                match_found <- TRUE
                match_field <- field
                match_value <- gene_info[[field]]
                break
              }
            }
            if (match_found) break
          }
        }
        
        if (!match_found) {
          next
        }
        
        cat("    ✓ Found matching gene in", match_field, ":", match_value, "\n")
        
        # Validate coordinates
        if (is.na(start_pos) || is.na(end_pos) || start_pos < 1 || end_pos < 1 ||
            start_pos > length(genome_seq) || end_pos > length(genome_seq)) {
          cat("    ERROR: Invalid coordinates - start:", start_pos, "end:", end_pos, "\n")
          next
        }
        
        # Extract sequence
        gene_seq <- Biostrings::getSeq(genome_seq, start = start_pos, end = end_pos)
        
        # Apply reverse complement if needed
        if (is_complement) {
          cat("    Applying reverse complement\n")
          gene_seq <- Biostrings::reverseComplement(gene_seq)
        }
        
        gene_seq_char <- as.character(gene_seq)
        
        # Create gene name for header
        gene_name <- match_value
        if (gene_name == "" || gene_name == "unknown") {
          gene_name <- paste0("feature_", i)
        }
        gene_name <- gsub("\\s+", "_", gene_name)
        
        # Create FASTA header
        strand_info <- ifelse(is_complement, "_complement", "")
        header <- paste0(">", genome_id, "_", gene_name, "_", feature_type, "_", 
                         start_pos, "..", end_pos, strand_info)
        
        # Add to results
        all_extracted_sequences <- c(all_extracted_sequences, header, gene_seq_char)
        
        cat("    ✓ Successfully extracted", gene_name, "sequence (", nchar(gene_seq_char), "bp)\n")
        cat("    Sequence preview:", substr(gene_seq_char, 1, 60), "...\n")
      }
    }
    
    if (length(all_extracted_sequences) == 0) {
      stop("No matching genes found in any specified feature types")
    }
    
    cat("    Total gene sequences extracted:", length(all_extracted_sequences) / 2, "\n")
    return(all_extracted_sequences)
    
  }, error = function(e) {
    stop("Gene extraction failed: ", e$message)
  })
}
# FIXED: Main function with proper GenBank file saving
get_gene_target_by_taxid <- function(taxid_file, max_per_taxid, output_dir, gene_synonyms, feature_types = c("CDS")) {
  
  if (missing(gene_synonyms) || length(gene_synonyms) == 0) {
    stop("gene_synonyms parameter is required and must contain at least one gene name")
  }
  
  if (!file.exists(taxid_file)) {
    stop("File ", taxid_file, " not found!")
  }
  
  taxids <- readLines(taxid_file)
  taxids <- taxids[taxids != ""]
  
  cat("Processing", length(taxids), "taxonomic IDs\n")
  cat("Getting up to", max_per_taxid, "sequences per taxid\n")
  cat("Searching for genes:", paste(gene_synonyms, collapse = ", "), "\n")
  cat("Looking for feature types:", paste(feature_types, collapse = ", "), "\n")
  
  # Initialize results tracking
  results_summary <- data.frame(
    taxid = character(),
    taxon_name = character(),
    sequences_found = integer(),
    sequences_retrieved = integer(),
    gene_target_extracted_from_genomes = integer(),
    stringsAsFactors = FALSE
  )
  
  # Create output directories
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  genbank_dir <- file.path(output_dir, "genbank_records")
  if (!dir.exists(genbank_dir)) {
    dir.create(genbank_dir)
  }
  
  # Process each taxonomy ID
  for (i in seq_along(taxids)) {
    taxid <- trimws(taxids[i])
    cat("Processing taxid", taxid, "(", i, "of", length(taxids), ")\n")
    
    tryCatch({
      
      search_query <- create_gene_search_query(taxid, gene_synonyms)
      cat("  Query:", search_query, "\n")
      
      # Search with retries
      search_result <- NULL
      attempt <- 1
      taxon_name <- "Unknown"
      
      while (is.null(search_result) && attempt <= 3) {
        tryCatch({
          search_result <- rentrez::entrez_search(
            db = "nuccore",
            term = search_query,
            retmax = max_per_taxid * 3
          )
          tax_summary <- rentrez::entrez_summary(db = "taxonomy", id = taxid)
          if (!is.null(tax_summary$scientificname)) {
            taxon_name <- tax_summary$scientificname
          }
          cat("  Taxon name:", taxon_name, "\n")
        }, error = function(e) {
          cat("  Search attempt", attempt, "failed:", e$message, "\n")
          Sys.sleep(2 * attempt)
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
      
      # Initialize variables
      final_sequences <- character()
      sequences_retrieved <- 0
      gene_target_extracted_count <- 0
      
      # Process sequences
      if (sequences_to_process > 0) {
        for (seq_id in search_result$ids) {
          if (sequences_retrieved >= max_per_taxid) break
          
          attempt_seq <- 1
          success_seq <- FALSE
          
          while (attempt_seq <= 3 && !success_seq) {
            tryCatch({
              seq_summary <- rentrez::entrez_summary(db = "nuccore", id = seq_id)
              seq_type <- determine_sequence_type(seq_summary$title)
              
              if (seq_type == "GENOMIC") {
                cat("    Processing genomic sequence:", seq_id, "\n")
                
                # Fetch GenBank record
                genbank_record <- rentrez::entrez_fetch(
                  db = "nuccore",
                  id = seq_id,
                  rettype = "gbwithparts",
                  retmode = "text"
                )
                
                # Validate genomic sequence contains target genes
                if (validate_genomic_gene_match(genbank_record, gene_synonyms)) {
                  cat("    ✓ Genomic sequence contains target genes\n")
                  
                  # FIXED: Save GenBank file permanently (not as temp file)
                  safe_taxon_name <- gsub("[^A-Za-z0-9]", "_", taxon_name)
                  gbk_filename <- file.path(genbank_dir, paste0(seq_id, "_", safe_taxon_name, ".gbk"))
                  writeLines(genbank_record, gbk_filename)
                  cat("    Saved GenBank file:", gbk_filename, "\n")
                  
                  # Extract gene sequences
                  gene_sequences <- extract_gene_target_with_geneviewer(gbk_filename, seq_id, gene_synonyms, feature_types)
                  
                  if (length(gene_sequences) > 0) {
                    final_sequences <- c(final_sequences, gene_sequences)
                    extracted_count <- length(gene_sequences) / 2
                    gene_target_extracted_count <- gene_target_extracted_count + extracted_count
                    sequences_retrieved <- sequences_retrieved + extracted_count
                    cat("      ✓ Extracted", extracted_count, "gene sequences from genome\n")
                  }
                  
                } else {
                  cat("    ✗ Genomic sequence does not contain target genes, skipping\n")
                }
                
              } else {
                # Regular gene sequence
                cat("    Validating", seq_type, "sequence:", seq_id, "\n")
                
                if (validate_gene_match(seq_summary, gene_synonyms)) {
                  cat("    ✓ Sequence matches target genes, fetching...\n")
                  
                  gene_sequence <- rentrez::entrez_fetch(
                    db = "nuccore",
                    id = seq_id,
                    rettype = "fasta"
                  )
                  
                  lines <- strsplit(gene_sequence, "\n")[[1]]
                  header <- lines[1]
                  sequence <- paste(lines[-1], collapse = "")
                  
                  new_header <- paste0(">", seq_id, "_taxid_", taxid, "_", seq_type, "_", gsub("^>", "", header))
                  
                  final_sequences <- c(final_sequences, new_header, sequence)
                  sequences_retrieved <- sequences_retrieved + 1
                  cat("    ✓ Retrieved", seq_type, "sequence\n")
                  
                } else {
                  cat("    ✗ Sequence does not match target genes, skipping\n")
                }
              }
              
              Sys.sleep(0.3)
              success_seq <- TRUE
              
            }, error = function(e) {
              cat("    Attempt", attempt_seq, "failed for sequence", seq_id, ":", e$message, "\n")
              Sys.sleep(0.5 * attempt_seq)
              attempt_seq <<- attempt_seq + 1
              if (attempt_seq > 3) {
                cat("    Failed to process sequence", seq_id, "after 3 attempts\n")
                write(paste("Failed to process sequence:", seq_id, "error:", e$message), 
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
        cat("  ✓ Retrieved", sequences_retrieved, "sequences total\n")
        cat("  Gene features extracted from genomes:", gene_target_extracted_count, "\n")
        cat("  Saved to:", output_file, "\n")
      } else {
        cat("  No sequences obtained\n")
      }
      
      # Add to summary
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
  fasta_files <- list.files(output_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    cat("No FASTA files found in", output_dir, "\n")
    return()
  }
  
  cat("Combining", length(fasta_files), "FASTA files into", combined_file, "\n")
  
  all_sequences <- character()
  
  for (file in fasta_files) {
    sequences <- readLines(file)
    all_sequences <- c(all_sequences, sequences)
  }
  
  writeLines(all_sequences, combined_file)
  
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

# Run the analysis
results <- get_gene_target_by_taxid(
  taxid_file = "macroscelidea_taxid.txt",
  max_per_taxid = 5,
  output_dir = "test_macroscelidea",
  gene_synonyms = coi_synonyms,
  feature_types = c("CDS")
)

# Combine all sequences
combine_sequences(
  output_dir = "test_macroscelidea",
  combined_file = "test_macroscelidea/gene_sequences.fasta"
)