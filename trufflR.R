#!/usr/bin/env RScript

# TrufflR: Extract specific gene sequences from NCBI's database for specific taxa
# Command line version

#---------------------------------------------------------------------#

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

#---------------------------------------------------------------------#
# Define command line options

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
  
  make_option(c("--nt-file"), type="character", default="combined_nucleotide_seqs.fna",
              help="Name for combined nucleotide file [default=%default]", 
              metavar="FILENAME"),
  
  make_option(c("--aa-file"), type="character", default="combined_aminoacid_seqs.faa",
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

# VALIDATION AND DIRECTORY HANDLING
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

# Handle output directory
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

#---------------------------------------------------------------------#
# Script subfunctions start here

# Extract accession from FASTA header
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

# Initialize search query
initialize_search_query <- function(taxid, gene_synonyms) {
  # This function creates the search query for NCBI and prepares gene patterns for filtering
  # It also retrieves taxon information for better organization
  
  # Get taxon name for better file organization and error handling
  taxon_name <- "Unknown"
  taxon_name_clean <- paste0("taxid_", taxid)  # fallback name
  
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
  })
  
  # Build search query combining taxid and gene synonyms
  gene_query_part <- paste(gene_synonyms, collapse = " OR ")
  full_query <- paste0("txid", taxid, "[Organism] AND (", gene_query_part, ")")
  
  # Define gene patterns for filtering (remove NCBI field tags and convert to lowercase)
  gene_patterns <- gsub("\\[.*?\\]", "", gene_synonyms)
  gene_patterns <- tolower(gene_patterns)
  
  return(list(
    taxon_name = taxon_name,
    taxon_name_clean = taxon_name_clean,
    full_query = full_query,
    gene_patterns = gene_patterns
  ))
}

# Initialize storage lists
initialize_storage_lists <- function() {
  # This function creates all the empty data structures needed to store results
  # Having this separate makes it clear what data we're collecting
  
  return(list(
    gb_str_dfs = list(),           # GenBank feature dataframes
    gb_seq_dfs = list(),           # FASTA sequence objects
    sequence_info = list(),        # Metadata about each sequence
    id_to_accession = list()       # Mapping from NCBI IDs to accessions
  ))
}

# Perform NCBI search
perform_ncbi_search <- function(full_query, retmax) {
  # This function handles the actual NCBI database search
  # Separating this makes it easier to modify search parameters or add error handling
  
  search_result <- rentrez::entrez_search(
    db = "nuccore",
    term = full_query,
    retmax = retmax
  )
  
  return(search_result$ids)
}

# Save raw GenBank and nucleotide FASTA sequences
save_raw_sequences <- function(ids, raw_files_folder, storage_lists) {
  # This function downloads and saves the raw GenBank and FASTA files
  # It also populates the storage lists with parsed data
  
  for (i in seq_along(ids)) {
    cat("Processing ID:", ids[i], "\n")
    
    # Fetch FASTA record first to get accession and title
    fasta_records <- rentrez::entrez_fetch(
      db = "nuccore",
      id = ids[i],
      rettype = "fasta",
      retmode = "text"
    )
    
    # Extract accession from FASTA header using existing helper function
    fasta_lines <- strsplit(fasta_records, "\n")[[1]]
    header_line <- fasta_lines[1]
    sequence_title <- gsub("^>", "", header_line)
    accession <- extract_accession(header_line)
    
    # Store the mapping from NCBI ID to accession for consistent naming
    storage_lists$id_to_accession[[ids[i]]] <- accession
    
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
    
    # Store sequence metadata with both ID and accession
    storage_lists$sequence_info[[ids[i]]] <- list(
      id = ids[i],
      accession = accession,
      title = sequence_title,
      type = determine_sequence_type(sequence_title)  # using existing helper function
    )
    
    # Parse GenBank file and store the dataframe
    tryCatch({
      gb_df <- geneviewer::read_gbk(filename_gb) %>% 
        geneviewer::gbk_features_to_df()
      storage_lists$gb_str_dfs[[ids[i]]] <- gb_df
      
      cat("GenBank features parsed for", accession, "\n")
      
    }, error = function(e) {
      cat("Error processing GenBank for ID", ids[i], ":", e$message, "\n")
    })
    
    # Load FASTA sequence into Biostrings object
    tryCatch({
      gb_fa <- Biostrings::readDNAStringSet(filename_fa)
      storage_lists$gb_seq_dfs[[ids[i]]] <- gb_fa
      
      cat("FASTA sequence loaded for", accession, "\n")
      
    }, error = function(e) {
      cat("Error processing FASTA for ID", ids[i], ":", e$message, "\n")
    })
  }
  
  return(storage_lists)
}

# Extract coordinates for target gene synonym genes
extract_gene_coordinates <- function(storage_lists, gene_patterns, feature_type) {
  # This function filters the GenBank features to find target genes
  # and extracts their coordinates for sequence extraction

  filtered_gb_str_dfs <- list()  # To store filtered dfs

  gene_coords <- lapply(names(storage_lists$gb_str_dfs), function(id) {
    df <- storage_lists$gb_str_dfs[[id]]
    if (!is.null(df) && nrow(df) > 0 && "gene" %in% colnames(df)) {
      # Filter by gene name using pattern matching
      filtered_df <- df[grepl(paste(gene_patterns, collapse = "|"),
                              tolower(df$gene), ignore.case = TRUE), ]
      # Optionally filter by feature type
      if (feature_type != "all" && "type" %in% colnames(filtered_df)) {
        filtered_df <- filtered_df[tolower(filtered_df$type) == tolower(feature_type), ]
      }
      # Save filtered df back to storage_lists
      filtered_gb_str_dfs[[id]] <- filtered_df
      # Return only the columns needed for coordinates
      available_cols <- intersect(c("gene", "strand", "start", "end", "type"), colnames(filtered_df))
      if (length(available_cols) > 0 && nrow(filtered_df) > 0) {
        return(filtered_df[, available_cols, drop = FALSE])
      }
    }
    # If no match, save empty df
    filtered_gb_str_dfs[[id]] <- data.frame()
    return(data.frame())
  })
  names(gene_coords) <- names(storage_lists$gb_str_dfs)
  # Remove empty dataframes
  gene_coords <- gene_coords[sapply(gene_coords, nrow) > 0]
  filtered_gb_str_dfs <- filtered_gb_str_dfs[sapply(filtered_gb_str_dfs, nrow) > 0]

  # Update storage_lists in the parent environment
  storage_lists$gb_str_dfs <- filtered_gb_str_dfs

  return(gene_coords)
}

# Save translated amino acid sequences for target genes
save_translated_sequences <- function(storage_lists, extracted_folder, gene_patterns) {
  # This function extracts protein translations from GenBank features and saves them as amino acid FASTA files
  
  all_translations <- list()
  
  for (id in names(storage_lists$gb_str_dfs)) {
    gb_df <- storage_lists$gb_str_dfs[[id]]
    accession <- storage_lists$id_to_accession[[id]]
    
    if (!is.null(gb_df) && "translation" %in% colnames(gb_df)) {
      # Filter for rows with non-empty translations
      translations_df <- gb_df[!is.na(gb_df$translation) & gb_df$translation != "", ]
      
      # Filter for target genes matching the search synonyms 
      if (nrow(translations_df) > 0 && "gene" %in% colnames(translations_df)) {
        translations_df <- translations_df[grepl(paste(gene_patterns, collapse = "|"), 
                                                tolower(translations_df$gene), ignore.case = TRUE), ]
      }
      
      if (nrow(translations_df) > 0) {
        # Create descriptive sequence names combining gene information
        seq_names <- paste0(
          accession, "|",
          ifelse(!is.na(translations_df$gene) & translations_df$gene != "", 
                 translations_df$gene, "unknown_gene"), "|",
          ifelse(!is.na(translations_df$type) & translations_df$type != "", 
                 translations_df$type, "unknown_type")
        )
        
        # Create AAStringSet from translations
        translation_seqs <- Biostrings::AAStringSet(translations_df$translation)
        names(translation_seqs) <- seq_names
        
        # Store in the combined list for return
        all_translations[[paste0("taxid_", id)]] <- translation_seqs
        
        # Save the translations as FASTA file
        protein_folder <- file.path(extracted_folder, "protein")
        if (!dir.exists(protein_folder)) dir.create(protein_folder, recursive = TRUE)
        translation_filename <- file.path(protein_folder, paste0(accession, "_aa.faa"))
        Biostrings::writeXStringSet(translation_seqs, translation_filename)
        cat("  Translations saved to:", translation_filename, "\n")
      }
    }
  }
  
  return(all_translations)
}

# Save subsectioned nucleotide sequences for target genes
save_gene_subsequences <- function(gene_coords, storage_lists, extracted_folder, taxon_name_clean, taxid, gene_patterns) {
  # This function extracts the actual gene sequences from the full genomic sequences based on the coordinates found in the GenBank features
  # The gene_coords should already be filtered, but we include gene_patterns for consistency
  
  for (id in names(gene_coords)) {
    if (id %in% names(storage_lists$gb_seq_dfs)) {
      coords <- gene_coords[[id]]
      dna_seq <- storage_lists$gb_seq_dfs[[id]]
      accession <- storage_lists$id_to_accession[[id]]
      
      cat("Extracting gene sequences for accession", accession, "(ID:", id, ")\n")
      
      # Create separate sequences for this accession
      id_extracted_sequences <- list()
      
      # ADDITIONAL SAFETY CHECK: Verify that the coordinates are for target genes
      # This ensures consistency even if gene_coords filtering had issues
      if ("gene" %in% colnames(coords)) {
        coords <- coords[grepl(paste(gene_patterns, collapse = "|"), 
                              tolower(coords$gene), ignore.case = TRUE), ]
      }
      
      # Skip if no matching genes remain after filtering
      if (nrow(coords) == 0) {
        cat("  No target genes found for accession", accession, "after filtering\n")
        next
      }
      
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
            
            # Handle reverse strand by taking reverse complement
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
        
        # Get the sequence type for this accession from stored metadata
        seq_type <- storage_lists$sequence_info[[id]]$type
        
        # Create clean FASTA headers: accession|gene_name|[sequence_type]
        clean_names <- character(length(id_extracted_sequences))
        for (k in seq_along(id_extracted_sequences)) {
          seq_name <- names(id_extracted_sequences)[k]
          # Extract gene name from the sequence name (format: accession_genename_featuretype_number)
          parts <- strsplit(seq_name, "_")[[1]]
          gene_name <- parts[2]  # Second part is the gene name
          
          # Create the new clean header format
          clean_names[k] <- paste0(accession, "|", gene_name, "|[", seq_type, "]")
        }
        names(gene_stringset) <- clean_names
        
        # Create filename with taxon name and accession
        output_filename <- paste0(accession, "_nt.fna")
        output_file <- file.path(extracted_folder, output_filename)
        
        nucleotide_folder <- file.path(extracted_folder, "nucleotide")
        if (!dir.exists(nucleotide_folder)) dir.create(nucleotide_folder, recursive = TRUE)
        output_file <- file.path(nucleotide_folder, output_filename)
        Biostrings::writeXStringSet(gene_stringset, output_file)
        cat("Gene sequences for accession", accession, "saved to:", output_file, "\n")
      } else {
        cat("No gene sequences could be extracted for accession", accession, "\n")
      }
    }
  }
}

# Save summary CSV and stdout information
save_summary_information <- function(storage_lists, taxid, taxon_name, taxon_name_clean, full_query, 
                                   feature_type, ids, gene_coords, output_folder, extracted_folder,
                                   overall_sequence_summary, gene_synonyms) {
  # This function creates comprehensive summary information about the analysis
  # including CSV files and text summaries
  # fix: was missing the gene_synonyms param
  
  # Create sequence type summary table for this taxid
  seq_summary_df <- do.call(rbind, lapply(storage_lists$sequence_info, function(x) {
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
  
  # Add to overall summary (this will be passed back to main function)
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
  
  # Print filtered GenBank dataframes for this taxid
  cat("Filtered GenBank features for taxid", taxid, ":\n")
  for (id in names(storage_lists$gb_str_dfs)) {
    gb_df <- storage_lists$gb_str_dfs[[id]]
    accession <- storage_lists$id_to_accession[[id]]
    
    if (!is.null(gb_df) && nrow(gb_df) > 0 && "gene" %in% colnames(gb_df)) {
      # Apply the same filtering as in coordinate extraction
      gene_patterns <- gsub("\\[.*?\\]", "", gene_synonyms)
      gene_patterns <- tolower(gene_patterns)
      
      filtered_gb_df <- gb_df[grepl(paste(gene_patterns, collapse = "|"), 
                                    tolower(gb_df$gene), ignore.case = TRUE), ]
      
      if (feature_type != "all" && "type" %in% colnames(filtered_gb_df)) {
        filtered_gb_df <- filtered_gb_df[tolower(filtered_gb_df$type) == tolower(feature_type), ]
      }
      
      cat("Filtered GenBank features for", accession, ":\n")
      if (nrow(filtered_gb_df) > 0) {
        print(filtered_gb_df)
      } else {
        cat("No matching features found for the specified genes and feature type.\n")
      }
      cat("\n")
    }
  }
  
  # Print FASTA sequences for this taxid
  cat("All FASTA sequences for taxid", taxid, ":\n")
  for (id in names(storage_lists$gb_seq_dfs)) {
    accession <- storage_lists$id_to_accession[[id]]
    cat("FASTA sequences for", accession, ":\n")
    print(storage_lists$gb_seq_dfs[[id]])
    cat("\n")
  }
  
  # Create analysis summary file
  summary_file <- file.path(output_folder, paste0(taxon_name_clean, "_taxid", taxid, "_analysis_summary.txt"))
  cat("Analysis Summary for Taxid:", taxid, "\n", file = summary_file)
  cat("Taxon name:", taxon_name, "\n", file = summary_file, append = TRUE)
  cat("Query:", full_query, "\n", file = summary_file, append = TRUE)
  cat("Feature type filter:", feature_type, "\n", file = summary_file, append = TRUE)
  cat("Total sequences processed:", length(ids), "\n", file = summary_file, append = TRUE)
  cat("GenBank files created:", length(storage_lists$gb_str_dfs), "\n", file = summary_file, append = TRUE)
  cat("FASTA files created:", length(storage_lists$gb_seq_dfs), "\n", file = summary_file, append = TRUE)
  cat("Gene matches found:", length(gene_coords), "\n", file = summary_file, append = TRUE)
  
  # Count total extracted sequences across all accessions
  total_extracted <- length(list.files(extracted_folder, pattern = paste0("*nt.fna")))
  cat("Separate FASTA files created:", total_extracted, "\n", file = summary_file, append = TRUE)
  
  # Add sequence type breakdown to summary
  cat("\nSequence Type Breakdown:\n", file = summary_file, append = TRUE)
  for (type in names(type_counts)) {
    cat("  ", type, ":", type_counts[type], "\n", file = summary_file, append = TRUE)
  }
  
  # Add accession mapping to summary
  cat("\nNCBI ID to Accession Mapping:\n", file = summary_file, append = TRUE)
  for (ncbi_id in names(storage_lists$id_to_accession)) {
    cat("  ", ncbi_id, "->", storage_lists$id_to_accession[[ncbi_id]], "\n", file = summary_file, append = TRUE)
  }
  
  return(list(
    seq_summary_df = seq_summary_df,
    overall_sequence_summary = overall_sequence_summary,
    summary_file = summary_file
  ))
}

#---------------------------------------------------------------------#

# Main function that orchestrates all the subfunctions
extract_gene_sequences <- function(taxid_file, gene_synonyms, feature_type = "all", retmax = 5, output_base_dir) {
  
  # Read taxids from file (keeping this in main function as it's setup logic)
  if (!file.exists(taxid_file)) {
    stop("Taxid file not found: ", taxid_file)
  }
  
  taxids <- readLines(taxid_file)
  taxids <- trimws(taxids)  # Remove whitespace
  taxids <- taxids[taxids != ""]  # Remove empty lines
  
  cat("Searching for", taxids, "\n")
  cat("Gene synonyms:", paste(gene_synonyms, collapse = ", "), "\n")
  cat("Feature type filter:", feature_type, "\n\n")
  
  # Initialize overall storage for all taxids
  all_gb_str_dfs <- list()
  all_gb_seq_dfs <- list()
  all_translations <- list()
  all_results <- list()
  overall_sequence_summary <- data.frame()
  
  # Process each taxid using the subfunctions
  for (taxid in taxids) {
    cat("=" , rep("=", 10), "=\n")
    cat("Processing Taxid:", taxid, "\n")
    cat("=" , rep("=", 10), "=\n")
    
    tryCatch({
      # 1. Initialize search query
      search_info <- initialize_search_query(taxid, gene_synonyms)
      
      cat("Taxon name:", search_info$taxon_name, "\n")
      
      # 3. Perform NCBI search
      ids <- perform_ncbi_search(search_info$full_query, retmax)
      
      if (length(ids) == 0) {
        cat("No sequences found for taxid", taxid, "\n\n")
        next
      }
      
      cat("Found", length(ids), "sequences for taxid", taxid, "\n")
      
      # Create output folder structure (keeping this here as it's setup logic)
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
      
      # 2. Initialize storage lists
      storage_lists <- initialize_storage_lists()
      
      # 4. Save raw GenBank and nucleotide FASTA sequences
      storage_lists <- save_raw_sequences(ids, raw_files_folder, storage_lists)
      
      # Add to overall storage
      for (id in names(storage_lists$gb_str_dfs)) {
        all_gb_str_dfs[[paste0(taxid, "_", id)]] <- storage_lists$gb_str_dfs[[id]]
      }
      for (id in names(storage_lists$gb_seq_dfs)) {
        all_gb_seq_dfs[[paste0(taxid, "_", id)]] <- storage_lists$gb_seq_dfs[[id]]
      }
      
      # 5. Extract coordinates for target gene synonym genes
      gene_coords <- extract_gene_coordinates(storage_lists, search_info$gene_patterns, feature_type)
      
      # 6. Save translated amino acid sequences for target genes
      taxid_translations <- save_translated_sequences(storage_lists, extracted_folder, search_info$gene_patterns) # fix: was passing gene_synonyms but should have passed gene_patterns
      all_translations <- c(all_translations, taxid_translations)
      
      # 7. Save subsectioned nucleotide sequences for target genes
      save_gene_subsequences(gene_coords, storage_lists, extracted_folder, 
                           search_info$taxon_name_clean, taxid, search_info$gene_patterns)
      
      # 8. Save summary CSV and stdout information
      summary_results <- save_summary_information(storage_lists, taxid, search_info$taxon_name, 
                                                 search_info$taxon_name_clean, search_info$full_query,
                                                 feature_type, ids, gene_coords, output_folder, 
                                                 extracted_folder, overall_sequence_summary, gene_synonyms)
      
      overall_sequence_summary <- summary_results$overall_sequence_summary
      
      # Store results for this taxid
      all_results[[taxid]] <- list(
        taxid = taxid,
        taxon_name = search_info$taxon_name,
        query = search_info$full_query,
        ids = ids,
        id_to_accession = storage_lists$id_to_accession,
        gb_str_dfs = storage_lists$gb_str_dfs,
        fasta_sequences = storage_lists$gb_seq_dfs,
        sequence_info = storage_lists$sequence_info,
        sequence_summary = summary_results$seq_summary_df,
        gene_coordinates = gene_coords,
        output_folder = output_folder,
        summary_file = summary_results$summary_file
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
  
  # Print overall summary (keeping this in main function as it's final reporting)
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
    
    # Save overall summary table
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

# Combine all nt sequences into one multifasta file
combine_nt_sequences <- function(extracted_folder, combined_file = "combined_seqs.fna") {
  # Find all *_nt.fna files recursively in extracted_folder
  fasta_files <- list.files(extracted_folder, pattern = "_nt\\.fna$", full.names = TRUE, recursive = TRUE)
  
  if (length(fasta_files) == 0) {
    cat("No FASTA files found in", extracted_folder, "\n")
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

# Combine all aa sequences into one multifasta file
combine_aa_sequences <- function(extracted_folder, combined_file = "combined_seqs.faa") {
  # Find all *_aa.faa files recursively in extracted_folder
  fasta_files <- list.files(extracted_folder, pattern = "_aa\\.faa$", full.names = TRUE, recursive = TRUE)
  
  if (length(fasta_files) == 0) {
    cat("No AA FASTA files found in", extracted_folder, "\n")
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

#---------------------------------------------------------------------#
# MAIN EXECUTION
# Update the function call to use the cleaned output_dir variable
cat("Starting TrufflR analysis...\n")

results <- extract_gene_sequences(
  taxid_file = opt$taxids,
  gene_synonyms = gene_synonyms, # Use the cleaned variable
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