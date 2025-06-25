library(httr)
library(jsonlite)
library(argparse)


parser <- ArgumentParser(description = "Convert Ensembl IDs to gene symbols via biotools.fr")

parser$add_argument("-i", "--input", required=TRUE, help="Path to input file containing Ensembl IDs")
parser$add_argument("-o", "--output", required=TRUE, help="Path to output TSV file")

args <- parser$parse_args()

input_file <- args$input
output_file <- args$output


read_ids_from_file <- function(file_path) {
  tryCatch({
    # Read the file, each line becomes an element in a character vector
    ids <- suppressWarnings(readLines(file_path))
    
    # Remove any leading/trailing whitespace
    ids <- trimws(ids)
    
    # Remove empty lines
    ids <- ids[ids != ""]
    
    return(ids)
  }, error = function(e) {
    # Handle errors gracefully
    message("Error reading file: ", e$message)
    return(NULL) # Or return an empty vector, or stop, depending on your needs.
  })
}

###
# Multiple IDs to convert - use a POST request
###
url = "https://biotools.fr/human/ensembl_symbol_converter/"
ids <- read_ids_from_file(input_file)

ids_json <- toJSON(ids)

body <- list(api=1, ids=ids_json)
r <- POST(url, body = body)

output = fromJSON(content(r, "text"), flatten=TRUE)

# Initialize an empty list to store the results
gene_info <- list()

# Loop through each ID in the output dictionary and create a data frame
for (id in names(output)) {
  gene_name <- output[[id]]
  
  # If gene_name is empty or missing, replace with the ID itself
  if (is.null(gene_name) || gene_name == "") {
    gene_name <- id
  }
  
  # Append the ID and gene name to the list
  gene_info <- append(gene_info, list(c(id, gene_name)))
}

df <- do.call(rbind, gene_info)


df <- as.data.frame(df)
colnames(df) <- c("Ensembl_ID", "Gene_Name")

# Add suffixes to duplicate gene names
df$Gene_Name <- ave(df$Gene_Name, df$Gene_Name, FUN = function(x) {
  if (length(x) == 1) {
    return(x)
  } else {
    return(paste0(x, "_", letters[seq_along(x)]))
  }
})


if (nrow(df) == 0) {
  message("No data to write. The output file will not be created.")
} else {
  write.table(df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}