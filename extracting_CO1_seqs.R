library(rentrez)

# Function to get up to N CO1 sequences per taxonomic order
get_co1_by_order <- function(taxid_file = "test.txt", max_per_order = 5) {
  
  # Read taxonomy IDs
  if (!file.exists(taxid_file)) {
    stop("File ", taxid_file, " not found!")
  }
  
  taxids <- readLines(taxid_file)
  taxids <- taxids[taxids != ""]  # Remove empty lines
  
  cat("Processing", length(taxids), "taxonomic orders\n")
  cat("Getting up to", max_per_order, "CO1 sequences per order\n\n")
  
  # Initialize results tracking
  results_summary <- data.frame(
    taxid = character(),
    order_name = character(),
    sequences_found = integer(),
    sequences_retrieved = integer(),
    stringsAsFactors = FALSE
  )
  
  # Create output directory
  output_dir <- "co1_by_order"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Process each taxonomy ID
  for (i in seq_along(taxids)) {
    taxid <- trimws(taxids[i])
    cat("Processing taxid", taxid, "(", i, "of", length(taxids), ")\n")
    
    tryCatch({
      # Search for CO1 sequences in this taxonomic group
      search_query <- paste0(
        "txid", taxid, "[Organism] AND ",
        "(CO1[Gene] OR COI[Gene] OR COX1[Gene] OR",
        "\"cytochrome c oxidase subunit 1\"[All Fields] OR ",
        "\"cytochrome oxidase subunit I\"[All Fields])"
      )
      
      cat("  Query:", search_query, "\n")
      
      # Search with limited results
      search_result <- entrez_search(
        db = "nuccore",
        term = search_query,
        retmax = max_per_order
      )
      
      sequences_found <- search_result$count
      sequences_to_get <- length(search_result$ids)
      
      cat("  Found:", sequences_found, "total sequences\n")
      cat("  Retrieving:", sequences_to_get, "sequences\n")
      
      # Get taxonomic name for this taxid
      order_name <- "Unknown"
      if (sequences_to_get > 0) {
        # Try to get taxonomic information
        tax_summary <- entrez_summary(db = "taxonomy", id = taxid)
        if (!is.null(tax_summary$scientificname)) {
          order_name <- tax_summary$scientificname
        }
        
        cat("  Order name:", order_name, "\n")
        
        # Fetch the sequences
        sequences <- entrez_fetch(
          db = "nuccore",
          id = search_result$ids,
          rettype = "fasta"
        )
        
        # Save to individual file
        output_file <- file.path(output_dir, paste0("taxid_", taxid, "_", 
                                                    gsub("[^A-Za-z0-9]", "_", order_name), 
                                                    ".fasta"))
        writeLines(sequences, output_file)
        
        cat("  Saved to:", output_file, "\n")
      } else {
        cat("  No sequences found\n")
      }
      
      # Add to summary
      results_summary <- rbind(results_summary, data.frame(
        taxid = taxid,
        order_name = order_name,
        sequences_found = sequences_found,
        sequences_retrieved = sequences_to_get,
        stringsAsFactors = FALSE
      ))
      
      cat("  ✓ Complete\n\n")
      
      # Be nice to NCBI servers
      Sys.sleep(0.5)
      
    }, error = function(e) {
      cat("  ✗ Error:", e$message, "\n\n")
      
      # Add error entry to summary
      results_summary <<- rbind(results_summary, data.frame(
        taxid = taxid,
        order_name = "ERROR",
        sequences_found = 0,
        sequences_retrieved = 0,
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
  cat("Average sequences per order:", round(mean(results_summary$sequences_retrieved), 2), "\n\n")
  
  # Show top orders by sequence count
  cat("Top orders by sequences found:\n")
  top_orders <- results_summary[order(-results_summary$sequences_found), ][1:min(10, nrow(results_summary)), ]
  print(top_orders[, c("order_name", "sequences_found", "sequences_retrieved")])
  
  cat("\nFiles saved in:", output_dir, "/\n")
  cat("Summary saved as: retrieval_summary.csv\n")
  
  return(results_summary)
}

# Function to combine all sequences into one file (optional)
combine_sequences <- function(output_dir = "co1_by_order", combined_file = "all_co1_sequences.fasta") {
  fasta_files <- list.files(output_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    cat("No FASTA files found in", output_dir, "\n")
    return()
  }
  
  cat("Combining", length(fasta_files), "FASTA files into", combined_file, "\n")
  
  # Read and combine all sequences
  all_sequences <- character()
  for (file in fasta_files) {
    sequences <- readLines(file)
    all_sequences <- c(all_sequences, sequences)
  }
  
  writeLines(all_sequences, combined_file)
  
  # Count total sequences
  total_seqs <- length(grep("^>", all_sequences))
  cat("Combined file contains", total_seqs, "sequences\n")
  cat("Saved as:", combined_file, "\n")
}

# Main execution
results <- get_co1_by_order(taxid_file = "test.txt", max_per_order = 5)
