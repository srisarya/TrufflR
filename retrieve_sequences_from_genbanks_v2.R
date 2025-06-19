# TrufflR: Extract specific gene sequences from NCBI's database for specific taxa
library(rentrez)     # For interfacing with NCBI Entrez databases
library(geneviewer)  # For parsing GenBank files
library(Biostrings)  # For sequence manipulation
library(seqinr)      # For sequence manipulation and reverse complement

# Function to extract GenBank accession from FASTA header
# @param fasta_header - The FASTA header line (starting with >)
# @return Character string with accession number or original ID if not found
extract_accession <- function(fasta_header) {
  # Remove the > symbol from header
  header_clean <- gsub("^>", "", fasta_header)
  
  # Look for GenBank accession pattern: 2 letters + 6 digits + . + 1 digit
  # This regex captures the accession at the beginning of the header
  accession_match <- regexpr("^[A-Z]{2}[0-9]{6}\\.[0-9]", header_clean)
  
  if (accession_match > 0) {
    # Extract the matched accession
    accession <- substr(header_clean, accession_match, 
                        accession_match + attr(accession_match, "match.length") - 1)
    return(accession)
  } else {
    # Fallback: try to get first part before space (common format)
    first_part <- strsplit(header_clean, " ")[[1]][1]
    return(first_part)
  }
}

# Determine sequence type based on title
# @param -  title - The sequence title/description
# @return - Character - string indicating sequence type

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

# Extract gene sequences from GenBank for multiple taxa
# @param - taxid_file - Path to text file containing taxids (one per line)
# @param - gene_synonyms - Vector of gene search terms with NCBI field tags
# @param - feature_type - Feature type to extract: "CDS", "gene", "rRNA", "tRNA", or "all" (default: "all")
# @param - retmax - Maximum number of sequences to retrieve per taxid (default: 5)
# @param - output_base_dir - Base directory for output folders (default: current directory)
# @return - df - containing results for each taxid

extract_gene_sequences <- function(taxid_file, gene_synonyms, feature_type = "all", retmax = 5, output_base_dir) {
  
  # Read taxids from file
  if (!file.exists(taxid_file)) {
    stop("Taxid file not found: ", taxid_file)
  }
  
  taxids <- readLines(taxid_file)
  taxids <- trimws(taxids)  # Remove whitespace
  taxids <- taxids[taxids != ""]  # Remove empty lines
  
  cat("Searching for", taxids, "\n")
  cat("Gene synonyms:", paste(gene_synonyms, collapse = ", "), "\n")
  cat("Feature type filter:", feature_type, "\n\n")
  
  all_gb_str_dfs <- list()
  all_gb_seq_dfs <- list()
  
  # Define gene patterns for filtering (remove NCBI field tags and convert to lowercase)
  gene_patterns <- gsub("\\[.*?\\]", "", gene_synonyms)
  gene_patterns <- tolower(gene_patterns)
  
  # Store results for all taxids
  all_results <- list()
  # Store sequence type summary across all taxids
  overall_sequence_summary <- data.frame()
  
  # Process each taxid
  for (taxid in taxids) {
    cat("=" , rep("=", 10), "=\n")
    cat("Processing Taxid:", taxid, "\n")
    cat("=" , rep("=", 10), "=\n")
    
    tryCatch({
      # Get taxon name
      taxon_name <- "Unknown"
      tryCatch({
        tax_summary <- rentrez::entrez_summary(db = "taxonomy", id = taxid)
        if (!is.null(tax_summary$scientificname)) {
          taxon_name <- tax_summary$scientificname
          # Clean taxon name for file naming (remove spaces and special characters)
          taxon_name_clean <- gsub("[^A-Za-z0-9_]", "_", taxon_name)
          taxon_name_clean <- gsub("_+", "_", taxon_name_clean)  # Replace multiple underscores with single
          taxon_name_clean <- gsub("^_|_$", "", taxon_name_clean)  # Remove leading/trailing underscores
        }
      }, error = function(e) {
        cat("Warning: Could not retrieve taxon name for taxid", taxid, ":", e$message, "\n")
        taxon_name_clean <- paste0("taxid_", taxid)
      })
      
      cat("Taxon name:", taxon_name, "\n")
      
      # Build search query
      gene_query_part <- paste(gene_synonyms, collapse = " OR ")
      full_query <- paste0("txid", taxid, "[Organism] AND (", gene_query_part, ")")
      
      # Search NCBI
      search_result <- rentrez::entrez_search(
        db = "nuccore",
        term = full_query,
        retmax = retmax
      )
      
      ids <- search_result$ids
      
      if (length(ids) == 0) {
        cat("No sequences found for taxid", taxid, "\n\n")
        next
      }
      
      cat("Found", length(ids), "sequences for taxid", taxid, "\n")
      
      # Create output folder structure
      output_folder <- file.path(output_base_dir)
      if (!dir.exists(output_folder)) {
        dir.create(output_folder, recursive = TRUE)
        cat("Created folder:", output_folder, "\n")
      }
      
      raw_files_folder <- file.path(output_folder, "raw_files")
      if (!dir.exists(raw_files_folder)) {
        dir.create(raw_files_folder)
      }
      
      extracted_folder <- file.path(output_folder)
      if (!dir.exists(extracted_folder)) {
        dir.create(extracted_folder)
      }
      
      # Initialize storage lists
      gb_str_dfs <- list()
      gb_seq_dfs <- list()
      sequence_info <- list()  # Store sequence metadata including titles and accessions
      id_to_accession <- list()  # Map NCBI IDs to accessions for consistent naming
      
      # Process each sequence ID
      for (i in seq_along(ids)) {
        cat("Processing ID:", ids[i], "\n")
        
        # Fetch FASTA record first to get accession
        fasta_records <- rentrez::entrez_fetch(
          db = "nuccore",
          id = ids[i],
          rettype = "fasta",
          retmode = "text"
        )
        
        # Extract accession from FASTA header
        fasta_lines <- strsplit(fasta_records, "\n")[[1]]
        header_line <- fasta_lines[1]
        sequence_title <- gsub("^>", "", header_line)
        accession <- extract_accession(header_line)
        
        # Store the mapping from NCBI ID to accession
        id_to_accession[[ids[i]]] <- accession
        
        cat("  Accession:", accession, "\n")
        
        # Fetch and save GenBank record using accession for filename
        genbank_records <- rentrez::entrez_fetch(
          db = "nuccore",
          id = ids[i],
          rettype = "gb",
          retmode = "text"
        )
        filename_gb <- file.path(raw_files_folder, paste0(accession, ".gb"))
        writeLines(genbank_records, filename_gb)
        
        # Save FASTA record using accession for filename
        filename_fa <- file.path(raw_files_folder, paste0(accession, ".fa"))
        writeLines(fasta_records, filename_fa)
        
        # Store sequence info with both ID and accession
        sequence_info[[ids[i]]] <- list(
          id = ids[i],
          accession = accession,
          title = sequence_title,
          type = determine_sequence_type(sequence_title)
        )
        
        # Parse GenBank file
        tryCatch({
          gb_df <- geneviewer::read_gbk(filename_gb) %>% 
            geneviewer::gbk_features_to_df()
          gb_str_dfs[[ids[i]]] <- gb_df
          all_gb_str_dfs[[paste0(taxid, "_", ids[i])]] <- gb_df
          
          # Print the GenBank dataframe
          cat("GenBank features for", accession, ":\n")
          print(gb_df)
          cat("\n")
          
        }, error = function(e) {
          cat("Error processing GenBank for ID", ids[i], ":", e$message, "\n")
        })
        
        # Load FASTA sequence
        tryCatch({
          gb_fa <- Biostrings::readDNAStringSet(filename_fa)
          gb_seq_dfs[[ids[i]]] <- gb_fa
          all_gb_seq_dfs[[paste0(taxid, "_", ids[i])]] <- gb_fa
          
          # Print the FASTA sequences
          cat("FASTA sequences for", accession, ":\n")
          print(gb_fa)
          cat("\n")
          
        }, error = function(e) {
          cat("Error processing FASTA for ID", ids[i], ":", e$message, "\n")
        })
      }
      
      # Print all GenBank dataframes for this taxid
      cat("All GenBank dataframes for taxid", taxid, ":\n")
      print(gb_str_dfs)
      cat("\n")
      
      # Print all FASTA sequences for this taxid  
      cat("All FASTA sequences for taxid", taxid, ":\n")
      print(gb_seq_dfs)
      cat("\n")
      
      # Create sequence type summary table for this taxid (now includes accessions)
      seq_summary_df <- do.call(rbind, lapply(sequence_info, function(x) {
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
      
      # Add to overall summary
      if (nrow(seq_summary_df) > 0) {
        overall_sequence_summary <- rbind(overall_sequence_summary, seq_summary_df)
      }
      
      # Print sequence type summary for this taxid
      cat("\nSequence Type Summary for Taxid", taxid, ":\n")
      type_counts <- table(seq_summary_df$Sequence_Type)
      for (type in names(type_counts)) {
        cat("  ", type, ":", type_counts[type], "\n")
      }
      cat("\n")
      
      # Filter dataframes for target genes and feature types
      gene_coords <- lapply(gb_str_dfs, function(df) {
        if (!is.null(df) && nrow(df) > 0 && "gene" %in% colnames(df)) {
          # Filter by gene name
          filtered_df <- df[grepl(paste(gene_patterns, collapse = "|"), 
                                  tolower(df$gene), ignore.case = TRUE), ]
          
          # Filter by feature type if specified
          if (feature_type != "all" && "type" %in% colnames(filtered_df)) {
            filtered_df <- filtered_df[tolower(filtered_df$type) == tolower(feature_type), ]
          }
          
          available_cols <- intersect(c("gene", "strand", "start", "end", "type"), colnames(filtered_df))
          if (length(available_cols) > 0) {
            return(filtered_df[, available_cols, drop = FALSE])
          }
        }
        return(data.frame())
      })
      
      # Remove empty dataframes
      gene_coords <- gene_coords[sapply(gene_coords, nrow) > 0]
      
      # Extract gene sequences and save separately for each ID (now using accessions)
      for (id in names(gene_coords)) {
        if (id %in% names(gb_seq_dfs)) {
          coords <- gene_coords[[id]]
          dna_seq <- gb_seq_dfs[[id]]
          accession <- id_to_accession[[id]]  # Get the accession for this ID
          
          cat("Extracting gene sequences for accession", accession, "(ID:", id, ")\n")
          
          # Create separate sequences for this accession
          id_extracted_sequences <- list()
          
          for (j in 1:nrow(coords)) {
            gene_name <- coords$gene[j]
            start_pos <- coords$start[j]
            end_pos <- coords$end[j]
            strand <- if("strand" %in% colnames(coords)) coords$strand[j] else "+"
            feature_info <- if("type" %in% colnames(coords)) coords$type[j] else "unknown"
            
            tryCatch({
              if (length(dna_seq) > 0) {
                full_seq <- dna_seq[[1]]
                gene_subseq <- Biostrings::subseq(full_seq, start = start_pos, end = end_pos)
                
                if (strand == "-" || strand == -1) {
                  gene_subseq <- Biostrings::reverseComplement(gene_subseq)
                }
                
                # Use accession in sequence name instead of NCBI ID
                seq_name <- paste0(accession, "_", gene_name, "_", feature_info, "_", j)
                id_extracted_sequences[[seq_name]] <- gene_subseq
                cat("  Extracted:", gene_name, "(", feature_info, ",", length(gene_subseq), "bp )\n")
              }
            }, error = function(e) {
              cat("  Error extracting sequence:", e$message, "\n")
            })
          }
          
          # Save extracted sequences for this accession if any were found
          if (length(id_extracted_sequences) > 0) {
            gene_stringset <- Biostrings::DNAStringSet(id_extracted_sequences)
            
            # Get the sequence type for this accession from our stored info
            seq_type <- sequence_info[[id]]$type
            
            # Create clean FASTA headers: accession|gene_name|[sequence_type]
            clean_names <- character(length(id_extracted_sequences))
            for (k in seq_along(id_extracted_sequences)) {
              seq_name <- names(id_extracted_sequences)[k]
              # Extract just the gene name from the original sequence name
              # Original format was: accession_genename_featuretype_number
              parts <- strsplit(seq_name, "_")[[1]]
              gene_name <- parts[2]  # Second part is the gene name
              
              # Create the new clean header format
              clean_names[k] <- paste0(accession, "|", gene_name, "|[", seq_type, "]")
            }
            names(gene_stringset) <- clean_names
            
            # Create filename with taxon name and accession (not NCBI ID)
            output_filename <- paste0(taxon_name_clean, "_taxid", taxid, "_", accession, "_extracted_genes.fasta")
            output_file <- file.path(extracted_folder, output_filename)
            
            Biostrings::writeXStringSet(gene_stringset, output_file)
            cat("Gene sequences for accession", accession, "saved to:", output_file, "\n")
          } else {
            cat("No gene sequences could be extracted for accession", accession, "\n")
          }
        }
      }
      
      # Create analysis summary file (now includes accession information)
      summary_file <- file.path(output_folder, paste0(taxon_name_clean, "_taxid", taxid, "_analysis_summary.txt"))
      cat("Analysis Summary for Taxid:", taxid, "\n", file = summary_file)
      cat("Taxon name:", taxon_name, "\n", file = summary_file, append = TRUE)
      cat("Query:", full_query, "\n", file = summary_file, append = TRUE)
      cat("Feature type filter:", feature_type, "\n", file = summary_file, append = TRUE)
      cat("Total sequences processed:", length(ids), "\n", file = summary_file, append = TRUE)
      cat("GenBank files created:", length(gb_str_dfs), "\n", file = summary_file, append = TRUE)
      cat("FASTA files created:", length(gb_seq_dfs), "\n", file = summary_file, append = TRUE)
      cat("Gene matches found:", length(gene_coords), "\n", file = summary_file, append = TRUE)
      
      # Count total extracted sequences across all accessions (pattern updated)
      total_extracted <- length(list.files(extracted_folder, pattern = paste0(taxon_name_clean, "_taxid", taxid, "_.*_extracted_genes.fasta")))
      cat("Separate FASTA files created:", total_extracted, "\n", file = summary_file, append = TRUE)
      
      # Add sequence type breakdown to summary
      cat("\nSequence Type Breakdown:\n", file = summary_file, append = TRUE)
      for (type in names(type_counts)) {
        cat("  ", type, ":", type_counts[type], "\n", file = summary_file, append = TRUE)
      }
      
      # Add accession mapping to summary
      cat("\nNCBI ID to Accession Mapping:\n", file = summary_file, append = TRUE)
      for (ncbi_id in names(id_to_accession)) {
        cat("  ", ncbi_id, "->", id_to_accession[[ncbi_id]], "\n", file = summary_file, append = TRUE)
      }
      
      # Store results for this taxid (now includes accession mapping)
      all_results[[taxid]] <- list(
        taxid = taxid,
        taxon_name = taxon_name,
        query = full_query,
        ids = ids,
        id_to_accession = id_to_accession,
        gb_str_dfs = gb_str_dfs,
        fasta_sequences = gb_seq_dfs,
        sequence_info = sequence_info,
        sequence_summary = seq_summary_df,
        gene_coordinates = gene_coords,
        output_folder = output_folder,
        summary_file = summary_file
      )
      
      cat("Completed processing taxid", taxid, "\n\n")
      
    }, error = function(e) {
      cat("Error processing taxid", taxid, ":", e$message, "\n\n")
      all_results[[taxid]] <- list(
        taxid = taxid,
        error = e$message
      )
    })
  }
  
  # Print overall summary
  cat("=" , rep("=", 15), "=\n")
  cat("OVERALL SUMMARY\n")
  cat("=" , rep("=", 15), "=\n")
  cat("Total taxids processed:", length(taxids), "\n")
  cat("Successful:", sum(sapply(all_results, function(x) is.null(x$error))), "\n")
  cat("Failed:", sum(sapply(all_results, function(x) !is.null(x$error))), "\n")
  
  # Print overall sequence type summary
  if (nrow(overall_sequence_summary) > 0) {
    cat("\nOVERALL SEQUENCE TYPE SUMMARY:\n")
    overall_type_counts <- table(overall_sequence_summary$Sequence_Type)
    for (type in names(overall_type_counts)) {
      cat("  ", type, ":", overall_type_counts[type], "\n")
    }
    
    # Save overall summary table (now includes accession column)
    return(list(
      results = all_results,
      gb_dataframes = all_gb_str_dfs,     # ALL dataframes across all taxids
      fasta_sequences = all_gb_seq_dfs,   # ALL sequences across all taxids  
      sequence_summary = overall_sequence_summary
    ))
    
    overall_summary_file <- file.path(output_base_dir, "overall_sequence_summary.csv")
    write.csv(overall_sequence_summary, overall_summary_file, row.names = FALSE)
    cat("\nOverall sequence summary saved to:", overall_summary_file, "\n")
    
    # Print detailed table
    cat("\nDETAILED SEQUENCE TABLE:\n")
    cat("=" , rep("=", 50), "=\n")
    print(overall_sequence_summary)
  }
  
  return(invisible(all_results))
}

# Combine all output sequences into one multifasta file
# @param - output_dir - directory for output multifasta (default: current directory)
# @param - combined_file - output multifasta (default: combined_seqs.fa)

combine_sequences <- function(output_dir, combined_file = "combined_seqs.fa") {
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

# Example usage:

# Define gene synonyms
coi_synonyms <- c("COI[Gene]",
                  "COX1[Gene]",
                  "cytochrome c oxidase subunit 1[Title]",
                  "cytochrome oxidase[All Fields]")

rRNA_16S_synonyms <- c("16S[Title]", 
                       "16S rRNA[Title]", 
                       "16S ribosomal RNA[Title]",
                       "small subunit ribosomal RNA[Title]",
                       "SSU rRNA[Title]")

rRNA_28S_synonyms <- c("28S[Title]", 
                       "28S rRNA[Title]", 
                       "28S ribosomal RNA[Title]",
                       "large subunit ribosomal RNA[Title]",
                       "LSU rRNA[Title]")

# Run analysis with feature type filtering
results <- extract_gene_sequences("taxids/proseriate_taxid.txt", coi_synonyms, feature_type = "CDS", retmax = 10, output_base_dir = "test/test_proseriate")
combine_sequences(output_dir = "tests/test_proseriate", combined_file = "test/test_proseriate/gene_sequences.fasta")

#main_analysis_result<-extract_gene_sequences("results/rhabditophora_family_results/", coi_synonyms, feature_type = "CDS", retmax = 10, output_base_dir = "results/rhabditophora_family_results/")
#combine_sequences(output_dir = "results/rhabditophora_family_results/", combined_file = "results/rhabditophora_family_results/gene_sequences.fasta")

#main_analysis_result<-extract_gene_sequences("metazoan_orders.txt", coi_synonyms, feature_type = "CDS", retmax = 3, output_base_dir = "results/metazoan_orders_results")
#combine_sequences(output_dir = "results/metazoan_orders_results", combined_file = "results/metazoan_orders_results/gene_sequences.fasta")

# results <- extract_gene_sequences("metazoans_5_taxids.txt", coi_synonyms, feature_type = "tRNA", retmax = 3, output_base_dir = "test_metazoan_tRNA")
# combine_sequences(output_dir = "test_metazoan_tRNA", combined_file = "test_metazoan_tRNA/gene_sequences.fasta")

# results_SSU <- extract_gene_sequences("metazoans_5_taxids.txt", rRNA_16S_synonyms, feature_type = "rRNA", retmax = 3, output_base_dir = "test_metazoan_SSU")
# combine_sequences(output_dir = "test_metazoan_SSU", combined_file = "test_metazoan_SSU/gene_sequences.fasta")

# results_LSU <- extract_gene_sequences("metazoans_5_taxids.txt", rRNA_28S_synonyms, feature_type = "rRNA", retmax = 3, output_base_dir = "test_metazoan_LSU")
# combine_sequences(output_dir = "test_metazoan_LSU", combined_file = "test_metazoan_LSU/gene_sequences.fasta")