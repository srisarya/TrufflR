#!/usr/bin/env Rscript

# TrufflR: Extract specific gene sequences from NCBI's database for specific taxa
# Command line version

# install required packages if not already installed

# CRAN packages
if (!require("optparse", quietly = TRUE))
    install.packages("optparse")
library(optparse)

if (!require("rentrez", quietly = TRUE))
    install.packages("rentrez")
library(rentrez)

if (!require("seqinr", quietly = TRUE))
    install.packages("seqinr")
library(seqinr)

if (!require("geneviewer", quietly = TRUE))
    install.packages("geneviewer")
library(geneviewer)

# Bioconductor manager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Bioconductor packages
if (!require("Biostrings", quietly = TRUE))
    BiocManager::install("Biostrings")
library(Biostrings)

# Define command line options
# Define command line options - CORRECTED VERSION
option_list <- list(
  make_option(c("-t", "--taxids"), type="character", 
              help="Path to text file containing taxids (one per line)", 
              metavar="FILE"),
  
  make_option(c("-g", "--genes"), type="character", 
              help="Comma-separated list of gene search terms (e.g., 'COI[Gene],COX1[Gene],cytochrome c oxidase subunit I[Gene]')", 
              metavar="GENES"),
  
  make_option(c("-f", "--feature-type"), type="character", default="all",
              help="Feature type to extract: CDS, gene, rRNA, tRNA, or all [default=%default]", 
              metavar="TYPE"),
  
  make_option(c("-r", "--retmax"), type="integer", default=5,
              help="Maximum number of sequences to retrieve per taxid [default=%default]", 
              metavar="NUMBER"),
  
  make_option(c("-o", "--output-dir"), type="character", default="trufflr_output",
              help="Output directory [default=%default]", 
              metavar="DIR"),
  
  make_option(c("-c", "--combine-nt"), action="store_true", default=FALSE,
              help="Combine all nucleotide sequences into one file"),
  
  make_option(c("-a", "--combine-aa"), action="store_true", default=FALSE,
              help="Combine all amino acid sequences into one file"),
  
  make_option(c("--nt-file"), type="character", default="combined_seqs.fna",
              help="Name for combined nucleotide file [default=%default]", 
              metavar="FILENAME"),
  
  make_option(c("--aa-file"), type="character", default="combined_seqs.faa",
              help="Name for combined amino acid file [default=%default]", 
              metavar="FILENAME"),
  
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output")
)
# Parse command line arguments
opt_parser <- OptionParser(
  option_list = option_list,
  description = "Extract specific gene sequences from NCBI's database for specific taxa",
  epilogue = paste(
    "Examples:",
    "  Rscript trufflr.R -t taxids.txt -g 'COI[Gene],COX1[Gene]' -o results",
    "  Rscript trufflr.R -t taxids.txt -g 'COI[Gene]' -f CDS -r 10 -c -a",
    "",
    "Gene search terms should include NCBI field tags like [Gene] or [All Fields]",
    sep = "\n"
  )
)

opt <- parse_args(opt_parser)

# IMPROVED VALIDATION AND DIRECTORY HANDLING
# Validate required arguments first
if (is.null(opt$taxids)) {
  cat("Error: Taxids file is required (use -t or --taxids)\n")
  print_help(opt_parser)
  quit(status = 1)
}

if (is.null(opt$genes)) {
  cat("Error: Gene search terms are required (use -g or --genes)\n")
  print_help(opt_parser)
  quit(status = 1)
}

# Handle output directory with improved logic
if (is.null(opt$`output-dir`)) {
    # This should not happen with default value, but just in case
    opt$`output-dir` <- "trufflr_output"
    cat("Warning: Output directory was NULL, using default: trufflr_output\n")
}

# Clean up the output directory path
output_dir <- trimws(as.character(opt$`output-dir`))

# Additional safety check
if (output_dir == "" || is.na(output_dir)) {
    output_dir <- "trufflr_output"
    cat("Warning: Output directory was empty, using default: trufflr_output\n")
}

# Create the directory if it doesn't exist
if (!dir.exists(output_dir)) {
    tryCatch({
        dir.create(output_dir, recursive = TRUE)
        if (opt$verbose) cat("Created output directory:", output_dir, "\n")
    }, error = function(e) {
        cat("Error: Could not create output directory '", output_dir, "': ", e$message, "\n")
        quit(status = 1)
    })
} else {
    if (opt$verbose) cat("Using existing output directory:", output_dir, "\n")
}

# Verify directory is writable
test_file <- file.path(output_dir, "test_write.tmp")
tryCatch({
    writeLines("test", test_file)
    unlink(test_file)  # Clean up test file
}, error = function(e) {
    cat("Error: Output directory '", output_dir, "' is not writable: ", e$message, "\n")
    quit(status = 1)
})

# Check if taxids file exists
if (!file.exists(opt$taxids)) {
  cat("Error: Taxids file not found:", opt$taxids, "\n")
  quit(status = 1)
}

# Parse gene synonyms from comma-separated string
gene_synonyms <- trimws(strsplit(opt$genes, ",")[[1]])

# Print configuration if verbose
if (opt$verbose) {
  cat("Configuration:\n")
  cat("  Taxids file:", opt$taxids, "\n")
  cat("  Gene synonyms:", paste(gene_synonyms, collapse = ", "), "\n")
  cat("  Feature type:", opt$`feature-type`, "\n")
  cat("  Max sequences per taxid:", opt$retmax, "\n")
  cat("  Output directory:", output_dir, "\n")
  cat("  Combine nucleotides:", opt$`combine-nt`, "\n")
  cat("  Combine amino acids:", opt$`combine-aa`, "\n")
  cat("\n")
}

# Function to extract GenBank accession from FASTA header
extract_accession <- function(fasta_header) {
  header_clean <- gsub("^>", "", fasta_header)
  accession_match <- regexpr("^[A-Z]{2}[0-9]{6}\\.[0-9]", header_clean)
  
  if (accession_match > 0) {
    accession <- substr(header_clean, accession_match, 
                        accession_match + attr(accession_match, "match.length") - 1)
    return(accession)
  } else {
    first_part <- strsplit(header_clean, " ")[[1]][1]
    return(first_part)
  }
}

# Determine sequence type based on title
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
  all_translations <- list()
  
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
          
          # Add this right after the translation extraction block:
          if (!is.null(gb_df) && "translation" %in% colnames(gb_df)) {
            # Filter for rows with non-empty translations
            translations_df <- gb_df[!is.na(gb_df$translation) & gb_df$translation != "", ]
            
            if (nrow(translations_df) > 0) {
              # Create sequence names combining gene info
              seq_names <- paste0(
                accession, "_",
                ifelse(!is.na(translations_df$gene) & translations_df$gene != "", 
                       translations_df$gene, "unknown_gene"), "_",
                ifelse(!is.na(translations_df$type) & translations_df$type != "", 
                       translations_df$type, "unknown_type"), "_",
                seq_len(nrow(translations_df))
              )
              
              # Create AAStringSet from translations
              translation_seqs <- Biostrings::AAStringSet(translations_df$translation)
              names(translation_seqs) <- seq_names
              
              # Store in the combined list
              all_translations[[paste0(taxid, "_", ids[i])]] <- translation_seqs
              
              # SAVE THE TRANSLATIONS AS FASTA FILE:
              translation_filename <- file.path(output_base_dir, paste0(accession, "_translations.fasta"))
              Biostrings::writeXStringSet(translation_seqs, translation_filename)
              cat("  Translations saved to:", translation_filename, "\n")
            }
          }
          
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
  
  # Print overall sequence type summary
  if (nrow(overall_sequence_summary) > 0) {
    cat("\nOVERALL SEQUENCE TYPE SUMMARY:\n")
    overall_type_counts <- table(overall_sequence_summary$Sequence_Type)
    for (type in names(overall_type_counts)) {
      cat("  ", type, ":", overall_type_counts[type], "\n")
    }
    
    # Save overall summary table (now includes accession column)
    overall_summary_file <- file.path(output_base_dir, "overall_sequence_summary.csv")
    write.csv(overall_sequence_summary, overall_summary_file, row.names = FALSE)
    cat("\nOverall sequence summary saved to:", overall_summary_file, "\n")
    
    # Print detailed table
    cat("\nDETAILED SEQUENCE TABLE:\n")
    cat("=" , rep("=", 50), "=\n")
    print(overall_sequence_summary)
  }
  
  return(list(
    results = all_results,
    gb_dataframes = all_gb_str_dfs,
    fasta_sequences = all_gb_seq_dfs,
    translations = all_translations,
    sequence_summary = overall_sequence_summary
  ))
}

# Combine all output sequences into one multifasta file
combine_nt_sequences <- function(output_dir, combined_file = "combined_seqs.fna") {
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

combine_aa_sequences <- function(output_dir, combined_file = "combined_seqs.faa") {
  fasta_files <- list.files(output_dir, pattern = "\\_translations.fasta$", full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    cat("No AA FASTA files found in", output_dir, "\n")
    return()
  }
  
  cat("Combining", length(fasta_files), "AA FASTA files into", combined_file, "\n")
  
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

# MAIN EXECUTION
# Update the function call to use the cleaned output_dir variable
cat("Starting TrufflR analysis...\n")

results <- extract_gene_sequences(
  taxid_file = opt$taxids,
  gene_synonyms = gene_synonyms,
  feature_type = opt$`feature-type`,  # Note: using backticks for hyphenated names
  retmax = opt$retmax,
  output_base_dir = output_dir  # Use the cleaned variable
)

# Combine sequences if requested
if (opt$`combine-nt`) {
  cat("\nCombining nucleotide sequences...\n")
  nt_output_file <- file.path(output_dir, opt$`nt-file`)
  combine_nt_sequences(output_dir, nt_output_file)
}

if (opt$`combine-aa`) {
  cat("\nCombining amino acid sequences...\n")
  aa_output_file <- file.path(output_dir, opt$`aa-file`)
  combine_aa_sequences(output_dir, aa_output_file)
}

cat("\nTrufflR analysis complete!\n")
cat("Results saved in:", output_dir, "\n")