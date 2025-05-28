# Enhanced CO1 Sequence Retrieval with Completeness Analysis
# This script retrieves CO1 sequences from NCBI for taxonomic orders and analyzes their completeness

# Load required libraries
library(rentrez)  # For interfacing with NCBI Entrez databases

# Function to get up to N CO1 sequences per taxonomic order with completeness analysis
get_co1_by_order <- function(taxid_file, max_per_order, output_dir) {
  
  # Read taxonomy IDs from input file
  if (!file.exists(taxid_file)) {
    stop("File ", taxid_file, " not found!")
  }
  
  # Read all lines from the taxonomy ID file
  taxids <- readLines(taxid_file)
  taxids <- taxids[taxids != ""]  # Remove empty lines using base R subsetting
  
  # Print initial information about the analysis
  cat("Processing", length(taxids), "taxonomic orders\n")
  cat("Getting up to", max_per_order, "CO1 sequences per order\n")
  cat("Analyzing sequence completeness...\n\n")
  
  # Initialize results tracking dataframe with completeness columns
  results_summary <- data.frame(
    taxid = character(),                    # NCBI taxonomy ID
    order_name = character(),               # Scientific name of the taxonomic order
    sequences_found = integer(),            # Total CO1 sequences found in database
    sequences_retrieved = integer(),        # Number of sequences actually downloaded
    complete_genomes = integer(),           # Count of complete genome sequences
    complete_cds = integer(),               # Count of complete CDS sequences
    incomplete_cds = integer(),             # Count of incomplete/partial CDS sequences
    other_sequences = integer(),            # Count of other sequence types
    stringsAsFactors = FALSE
  )
  
  # Create output directory for storing FASTA files
  output_dir <- output_dir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Process each taxonomy ID sequentially
  for (i in seq_along(taxids)) {
    taxid <- trimws(taxids[i])  # Remove leading/trailing whitespace
    cat("Processing taxid", taxid, "(", i, "of", length(taxids), ")\n")
    
    # Use tryCatch for error handling - continue processing other taxids if one fails
    tryCatch({
      
      # Construct base search query for CO1 sequences
      # This searches for various CO1 gene synonyms and descriptions
      base_co1_query <- paste0(
        "txid", taxid, "[Organism] AND ",                    # Limit to specific taxonomic group
        "(CO1[Gene] OR COI[Gene] OR COX1[Gene] OR",          # Gene name variations
        "\"cytochrome c oxidase subunit 1\"[All Fields] OR ", # Full protein name
        "\"cytochrome oxidase subunit I\"[All Fields])"       # Alternative protein name
      )
      
      cat("  Base query:", base_co1_query, "\n")
      
      # Search for different types of sequences to analyze completeness
      
      # 1. Search for complete genomes containing CO1
      complete_genome_query <- paste0(base_co1_query, " AND complete genome[Title]")
      complete_genomes_result <- rentrez::entrez_search(
        db = "nuccore",                    # NCBI nucleotide database
        term = complete_genome_query,
        retmax = 0                         # Only get count, not actual IDs
      )
      complete_genomes_count <- complete_genomes_result$count
      
      # 2. Search for complete CDS sequences
      complete_cds_query <- paste0(base_co1_query, " AND (complete cds[Title] OR complete CDS[Title])")
      complete_cds_result <- rentrez::entrez_search(
        db = "nuccore",
        term = complete_cds_query,
        retmax = 0                         # Only get count, not actual IDs
      )
      complete_cds_count <- complete_cds_result$count
      
      # 3. Search for incomplete/partial CDS sequences
      incomplete_cds_query <- paste0(base_co1_query, " AND (partial[Title] OR incomplete[Title] OR partial cds[Title] OR partial CDS[Title])")
      incomplete_cds_result <- rentrez::entrez_search(
        db = "nuccore",
        term = incomplete_cds_query,
        retmax = 0                         # Only get count, not actual IDs
      )
      incomplete_cds_count <- incomplete_cds_result$count
      
      # 4. Main search to get actual sequence IDs for download
      search_result <- rentrez::entrez_search(
        db = "nuccore",
        term = base_co1_query,
        retmax = max_per_order             # Limit number of sequences to retrieve
      )
      
      # Calculate sequence counts
      sequences_found <- search_result$count                    # Total available sequences
      sequences_to_get <- length(search_result$ids)       # Number we'll actually download
      # Calculate "other" sequences (those not classified as complete/incomplete)
      other_sequences_count <- max(0, sequences_found - complete_genomes_count - complete_cds_count - incomplete_cds_count)
      
      # Print summary of what was found for this taxonomic group
      cat("  Found:", sequences_found, "total CO1 sequences\n")
      cat("    - Complete genomes:", complete_genomes_count, "\n")
      cat("    - Complete CDS:", complete_cds_count, "\n")
      cat("    - Incomplete CDS:", incomplete_cds_count, "\n")
      cat("    - Other sequences:", other_sequences_count, "\n")
      cat("  Retrieving:", sequences_to_get, "sequences\n")
      
      # Get taxonomic name for this taxid and download sequences if any found
      order_name <- "Unknown"
      if (sequences_to_get > 0) {
        
        # Retrieve taxonomic information from NCBI taxonomy database
        tax_summary <- rentrez::entrez_summary(db = "taxonomy", id = taxid)
        if (!is.null(tax_summary$scientificname)) {
          order_name <- tax_summary$scientificname
        }
        
        cat("  Order name:", order_name, "\n")
        
        # Download the actual FASTA sequences
        sequences <- rentrez::entrez_fetch(
          db = "nuccore",
          id = search_result$ids,
          rettype = "fasta"                # Get sequences in FASTA format
        )
        
        # Create safe filename by replacing non-alphanumeric characters with underscores
        safe_order_name <- gsub("[^A-Za-z0-9]", "_", order_name)
        output_file <- file.path(output_dir, paste0("taxid_", taxid, "_", safe_order_name, ".fasta"))
        
        # Write sequences to file
        writeLines(sequences, output_file)
        
        cat("  Saved to:", output_file, "\n")
      } else {
        cat("  No sequences found\n")
      }
      
      # Add results for this taxonomic group to summary dataframe
      results_summary <- rbind(results_summary, data.frame(
        taxid = taxid,
        order_name = order_name,
        sequences_found = sequences_found,
        sequences_retrieved = sequences_to_get,
        complete_genomes = complete_genomes_count,
        complete_cds = complete_cds_count,
        incomplete_cds = incomplete_cds_count,
        other_sequences = other_sequences_count,
        stringsAsFactors = FALSE
      ))
      
      cat("  ✓ Complete\n\n")
      
      # Be respectful to NCBI servers - wait between requests
      # Longer delay due to multiple queries per taxonomic ID
      Sys.sleep(1.0)
      
    }, error = function(e) {
      # Handle errors gracefully - log error but continue with next taxonomic ID
      cat("  ✗ Error:", e$message, "\n\n")
      
      # Add error entry to summary with zeros for all counts
      results_summary <<- rbind(results_summary, data.frame(
        taxid = taxid,
        order_name = "ERROR",
        sequences_found = 0,
        sequences_retrieved = 0,
        complete_genomes = 0,
        complete_cds = 0,
        incomplete_cds = 0,
        other_sequences = 0,
        stringsAsFactors = FALSE
      ))
    })
  }
  
  # Save comprehensive summary to CSV file
  write.csv(results_summary, file.path(output_dir, "retrieval_summary.csv"), row.names = FALSE)
  
  # Show top performing orders with completeness breakdown
  cat("Top orders by total sequences found:\n")
  # Sort by sequences_found in descending order and take top 10
  top_orders <- results_summary[order(-results_summary$sequences_found), ][1:min(10, nrow(results_summary)), ]
  # Print selected columns showing completeness breakdown
  print(top_orders[, c("order_name", "sequences_found", "complete_genomes", "complete_cds", "incomplete_cds", "other_sequences")])
  
  cat("\nFiles saved in:", output_dir, "/\n")
  cat("Enhanced summary saved as: retrieval_summary.csv\n")
  
  # Return the results dataframe for further analysis
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

# Main execution - run the analysis
results <- get_co1_by_order(taxid_file = "proseriate_taxid.txt", max_per_order = 10, output_dir = "proseriate_co1")
combine_sequences("proseriate_co1", "proseriate_co1s.fa")
