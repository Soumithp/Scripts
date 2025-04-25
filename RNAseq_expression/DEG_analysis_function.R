
# Function to perform DEG analysis with two-factor design and HTML output
runDEGanalysis <- function(
    gct_file,           # Input GCT file with raw counts
    metadata_file,      # Metadata text file with conditions
    output_html = "deg_results.html"  # Output HTML file
) {
  # Step 1: Read the GCT file
  d <- read.delim(gct_file, header = TRUE, skip = 2, check.names = FALSE)
  
  # Extract gene_id and gene_name
  gene_ids <- d$gene_id
  gene_names <- d$gene_name
  
  # Drop gene_id and gene_name columns from count matrix
  count_matrix <- d[, !(names(d) %in% c("gene_id", "gene_name"))]
  
  # Ensure counts are numeric and handle non-integers
  count_matrix <- as.matrix(count_matrix)
  if (!all(count_matrix == floor(count_matrix), na.rm = TRUE)) {
    warning("Some counts are not integers. Rounding to nearest integer.")
    count_matrix <- round(count_matrix)
  }
  
  # Filter out genes with all zeros or all NA
  keep <- rowSums(count_matrix > 0, na.rm = TRUE) > 0  # At least one non-zero value
  count_matrix <- count_matrix[keep, ]
  gene_ids <- gene_ids[keep]
  gene_names <- gene_names[keep]
  
  # Report how many genes were filtered
  cat("Total genes in GCT:", length(d$gene_id), "\n")
  cat("Genes retained after filtering (non-zero, non-NA):", length(gene_ids), "\n")
  
  # Set row names to gene_ids (DESeq2 requires unique row names)
  if (any(duplicated(gene_ids))) {
    warning("Duplicate gene_ids found. Making them unique.")
    gene_ids_unique <- make.unique(as.character(gene_ids))
  } else {
    gene_ids_unique <- gene_ids
  }
  rownames(count_matrix) <- gene_ids_unique
  
  # Step 2: Read the metadata file
  col_data <- read.delim(metadata_file, header = TRUE, row.names = 1)
  if (!all(c("sex", "treatment") %in% colnames(col_data))) {
    stop("Metadata file must contain 'sex' and 'treatment' columns.")
  }
  
  # Convert to factors
  col_data$sex <- factor(col_data$sex)
  col_data$treatment <- factor(col_data$treatment)
  
  # Ensure sample names match
  if (!all(colnames(count_matrix) %in% rownames(col_data))) {
    stop("Sample names in GCT file do not match metadata row names.")
  }
  col_data <- col_data[match(colnames(count_matrix), rownames(col_data)), ]
  
  # Step 3: Create DESeqDataSet object with two-factor design
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = col_data,
    design = ~ sex + treatment + sex:treatment
  )
  
  # Step 4: Run DESeq2 analysis
  dds <- DESeq(dds)
  
  # Print resultsNames for debugging
  cat("Available resultsNames:\n")
  print(resultsNames(dds))
  
  # Step 5: Extract results for default contrasts
  results_list <- list()
  treatments <- levels(col_data$treatment)
  ref_treatment <- treatments[1]  # Assume first level (e.g., "CDA") is reference
  
  # Treatment comparisons within Female
  for (i in 1:(length(treatments)-1)) {
    for (j in (i+1):length(treatments)) {
      contrast_name <- paste("F_", treatments[j], "_vs_", treatments[i], sep="")
      results_list[[contrast_name]] <- results(dds, contrast = c("treatment", treatments[j], treatments[i]))
    }
  }
  
  # Sex comparisons within each treatment
  for (treat in treatments) {
    contrast_name <- paste("M_vs_F_", treat, sep="")
    if (treat == ref_treatment) {
      results_list[[contrast_name]] <- results(dds, contrast = c("sex", "M", "F"))
    } else {
      interaction_term <- paste("sexM.treatment", treat, sep="")
      if (interaction_term %in% resultsNames(dds)) {
        results_list[[contrast_name]] <- results(dds, contrast = list("sex_M_vs_F", interaction_term))
      } else {
        warning(paste("Interaction term", interaction_term, "not found in resultsNames. Skipping."))
      }
    }
  }
  
  # Step 6: Generate HTML output
  html_content <- "<!DOCTYPE html>\n<html>\n<head>\n<title>DEG Results</title>\n"
  html_content <- paste0(html_content, "<style>\n",
                         "table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }\n",
                         "th, td { border: 1px solid black; padding: 8px; text-align: left; }\n",
                         "th { background-color: #f2f2f2; }\n",
                         "</style>\n</head>\n<body>\n<h1>DEG Analysis Results</h1>\n")
  
  for (name in names(results_list)) {
    res <- results_list[[name]]
    res_df <- as.data.frame(res)
    res_df <- cbind(gene_id = gene_ids, gene_name = gene_names, res_df)
    res_df <- res_df[order(res_df$padj), ]
    
    # Round numeric columns for better readability
    res_df$baseMean <- round(res_df$baseMean, 2)
    res_df$log2FoldChange <- round(res_df$log2FoldChange, 3)
    res_df$lfcSE <- round(res_df$lfcSE, 3)
    res_df$stat <- round(res_df$stat, 3)
    res_df$pvalue <- format.pval(res_df$pvalue, digits = 3)
    res_df$padj <- format.pval(res_df$padj, digits = 3)
    
    # Add table to HTML content
    html_content <- paste0(html_content, "<h2>", name, "</h2>\n")
    html_content <- paste0(html_content, kable(res_df, format = "html", row.names = FALSE) %>% 
                             kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE))
  }
  
  html_content <- paste0(html_content, "</body>\n</html>")
  
  # Write to HTML file
  writeLines(html_content, output_html)
  cat("DEG results written to:", output_html, "\n")
  
  # Return the DESeqDataSet object
  return(dds)
}
