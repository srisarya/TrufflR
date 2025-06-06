# TrufflR: Extract specific gene sequences from NCBI's database for specific taxa

# 1. create_gene_search_query() - Constructs NCBI search queries by combining
#    user-provided gene search terms with taxonomic constraints

# 2. extract_feature_location() - Helper function that parses GenBank feature 
#    location strings to extract start/end coordinates and strand information

# 3. extract_gene_target_from_genbank() - Core function that parses complete 
#    GenBank records, finds CDS features matching target genes, extracts their
#    sequences from the genome, and applies reverse complement if needed

# 4. get_gene_target_by_order() - Main orchestrating function that:
#    - Reads taxonomic IDs from input file
#    - For each taxon: searches NCBI, downloads sequences/genomes
#    - Calls extraction functions to get target genes
#    - Saves individual FASTA files and creates summary statistics

# 5. combine_sequences() - Final utility that merges all individual FASTA files
#    into one master file for downstream analysis

# USER INPUTS:
# - taxid_file: Text file containing NCBI taxonomic IDs (one per line)
# - gene_synonyms: Vector of search terms with NCBI field tags (e.g., "COI[Gene]")
# - max_per_order: Maximum sequences to retrieve per taxonomic group
# - output_dir: Directory path for saving results

# OUTPUTS:
# - Individual FASTA files per taxonomic group
# - Combined master FASTA file with all sequences
# - CSV summary with retrieval statistics
# - GenBank records folder with complete genome files
# - Error log for troubleshooting failed retrievals


# Load required libraries
library(rentrez)     # For interfacing with NCBI Entrez databases
library(Biostrings)  # For sequence manipulation

# Functions
create_gene_search_query <- function(taxid, gene_synonyms) {
  # Use gene synonyms exactly as provided by user (no additional field tags)
  # Combine all conditions with OR
  gene_query_part <- paste(gene_synonyms, collapse = " OR ")
  
  # Complete query - only add the taxid part
  full_query <- paste0("txid", taxid, "[Organism] AND (", gene_query_part, ")")
  
  return(full_query)
}

# Keep the original function as fallback
extract_gene_target_from_genbank <- function(genbank_text, genome_id, gene_synonyms) {
  # Combine GenBank text into a single string for easier handling
  gb_str <- paste(genbank_text, collapse = "\n")
  
  # Extract raw sequence from ORIGIN section
  origin_match <- regexpr("ORIGIN[^\n]*\n([a-z0-9 \n]*)//", gb_str, perl = TRUE)
  if (origin_match == -1) stop("No ORIGIN section found")
  
  origin_seq_raw <- regmatches(gb_str, origin_match)
  origin_lines <- unlist(strsplit(origin_seq_raw, "\n"))
  origin_seq <- tolower(gsub("[^acgt]", "", paste(origin_lines, collapse = "")))
  
  # Split text into lines and locate CDS features
  lines <- unlist(strsplit(gb_str, "\n"))
  cds_indices <- grep("^\\s{5}CDS\\s+", lines)
  
  for (i in cds_indices) {
    location_line <- lines[i]
    location_str <- sub("^\\s{5}CDS\\s+", "", location_line)
    
    # Capture qualifiers until next feature (starts with 5 spaces but not continuation)
    qualifiers <- c()
    j <- i + 1
    while (j <= length(lines) && (grepl("^\\s{21}", lines[j]) || grepl("^\\s{5}/", lines[j]))) {
      qualifiers <- c(qualifiers, lines[j])
      j <- j + 1
    }
    
    # Collapse wrapped lines (multiline qualifiers)
    qualifier_text <- paste(qualifiers, collapse = "")
    # Extract all gene and product qualifiers
    gene_matches <- regmatches(qualifier_text, gregexpr("/(gene|product)=['\"][^'\"]+['\"]", qualifier_text, perl = TRUE))[[1]]
    gene_values <- gsub("^/(gene|product)=['\"]|['\"]$", "", gene_matches)
    
    # Match any synonym
    if (any(tolower(gene_values) %in% tolower(gene_synonyms))) {
      # Handle complement() and join()
      is_complement <- grepl("^complement\\(", location_str)
      location_clean <- gsub("[^0-9\\.]", "", location_str)
      coords <- as.numeric(unlist(strsplit(location_clean, "\\.\\.")))
      
      if (length(coords) != 2) next  # Skip malformed location
      
      seq_sub <- substr(origin_seq, coords[1], coords[2])
      if (is_complement) {
        seq_sub <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq_sub)))
      }
      
      header <- paste0(">gene_target_from_", genome_id, "_", coords[1], "..", coords[2])
      return(c(header, seq_sub))
    }
  }
  
  stop("No matching gene found.")
}

# Renamed and refactored main function
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
              
              # Check if this is a complete genome
              if (!is.null(seq_summary$title) && grepl("complete genome", seq_summary$title, ignore.case = TRUE)) {
                
                # This is a complete genome - extract gene features using geneviewer
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
                
                # Extract gene sequences using geneviewer (with fallback)
                gene_sequences <- extract_gene_target_with_geneviewer(genbank_file, seq_id, gene_synonyms, feature_types)
                
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

# Example for rRNA searches
rrna_16s_synonyms <- c(
  "16S[Title]", 
  "16S rRNA[Title]", 
  "16S ribosomal RNA[rRNA]",
  "small subunit ribosomal RNA[Title]",
  "SSU rRNA[Title]"
)

# Run the analysis with COI synonyms (note the renamed function and parameter)
results <- get_gene_target_by_taxid(
  taxid_file = "proseriate_taxid.txt",
  max_per_taxid = 2,  # Renamed parameter
  output_dir = "test_proseriate_CO1",
  gene_synonyms = coi_synonyms,
  feature_types = c("CDS")  # Default for protein-coding genes
)

# Example for rRNA extraction
# results_rrna <- get_gene_target_by_taxid(
#   taxid_file = "metazoans_5_taxids.txt",
#   max_per_taxid = 2,
#   output_dir = "metazoans_5_taxids_16S",
#   gene_synonyms = rrna_16s_synonyms,
#   feature_types = c("rRNA")  # Specify rRNA features
# )

# Example for mixed feature types
# results_mixed <- get_gene_target_by_taxid(
#   taxid_file = "metazoans_5_taxids.txt",
#   max_per_taxid = 2,
#   output_dir = "metazoans_5_taxids_mixed",
#   gene_synonyms = c("COI", "16S"),
#   feature_types = c("CDS", "rRNA")  # Search both feature types
# )

# Combine all sequences into one file
combine_sequences(
  output_dir = "test_proseriate_CO1",
  combined_file = "test_proseriate_CO1/gene_sequences.fasta"
)