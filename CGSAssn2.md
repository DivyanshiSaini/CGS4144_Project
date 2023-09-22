#load library
BiocManager::install("org.Bt.eg.db") # Install cow genome annotations database

#lib
library("org.Bt.eg.db") # Load cow genome annotations

# We will need this so we can use the pipe: %>%
library(magrittr)


# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}


# Define the file path to the data directory
# Replace with the path of the folder the files will be in
data_dir <- file.path("AnalysisFolderNew","data", "SRP159810")

# Declare the file path to the gene expression matrix file
# inside directory saved as `data_dir`
# Replace with the path to your dataset file
data_file <- file.path(data_dir, "SRP159810.tsv")

# Declare the file path to the metadata file
# inside the directory saved as `data_dir`
# Replace with the path to your metadata file
metadata_file <- file.path(data_dir, "metadata_SRP159810.tsv")


# Check if the gene expression matrix file is at the path stored in `data_file`
file.exists(data_file)

# Check if the metadata file is at the file path stored in `metadata_file`
file.exists(metadata_file)

# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)


# Read in data TSV file
expression_df <- readr::read_tsv(data_file) %>%

# Tuck away the Gene ID column as row names
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)


# Bring back the "Gene" column in preparation for mapping
expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")

# Map Ensembl IDs to their associated Entrez IDs
mapped_list <- mapIds(
  org.Bt.eg.db, # Replace with annotation package for your organism
  keys = expression_df$Gene,
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "SYMBOL", # The type of gene identifiers you would like to map to
  multiVals = "list"
)

# Let's use the `head()` function for a preview of our mapped list
head(mapped_list)




# Let's make our list a bit more manageable by turning it into a data frame
mapped_df <- mapped_list %>%
  tibble::enframe(name = "Ensembl", value = "Symbol") %>%
  
# enframe() makes a `list` column; we will simplify it with unnest()
# This will result in one row of our data frame per list item
  tidyr::unnest(cols = Symbol)
  
  head(mapped_df)
  
  
# `maxsum = 10` limits the summary to 10 distinct values
summary(as.factor(mapped_df$Symbol), maxsum = 10)



multi_mapped <- mapped_df %>%
# Let's count the number of times each Ensembl ID appears in `Ensembl` column
  dplyr::count(Ensembl, name = "entrez_id_count") %>%
# Arrange by the genes with the highest number of Entrez IDs mapped
  dplyr::arrange(desc(entrez_id_count))

# Let's look at the first 6 rows of our `multi_mapped` object
head(multi_mapped)



collapsed_mapped_df <- mapped_df %>%
# Group by Ensembl IDs
  dplyr::group_by(Ensembl) %>%
# Collapse the Symbol IDs `mapped_df` into one column named `all_symbol_ids`
  dplyr::summarize(all_symbol_ids = paste(Symbol, collapse = ";"))


collapsed_mapped_df %>%
# Filter `collapsed_mapped_df` to include only the rows where
# `all_symbol_ids` values include the ";" character --
# these are the rows with multiple mapped values
  dplyr::filter(stringr::str_detect(all_symbol_ids, ";")) %>%
# We only need a preview here
  head()



final_mapped_df <- data.frame(
  "first_mapped_symbol_id" = mapIds(
    org.Bt.eg.db, # Replace with annotation package for your organism
    keys = expression_df$Gene,
    keytype = "ENSEMBL", # Replace with the gene identifiers used in your data
    column = "SYMBOL", # The type of gene identifiers you would like to map to
    multiVals = "first" # Keep only the first mapped value for each Ensembl ID
  )
) %>%
# Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
# Add the multiple mappings data from `collapsed_mapped_df` using Ensembl IDs
  dplyr::inner_join(collapsed_mapped_df, by = "Ensembl") %>%
# Now let's add on the rest of the expression data
  dplyr::inner_join(expression_df, by = c("Ensembl" = "Gene"))


#next
final_mapped_df %>%
  # Filter `final_mapped_df` to rows with multiple mapped values
  dplyr::filter(stringr::str_detect(all_symbol_ids, ";")) %>%
  head()


# Write mapped and annotated data frame to output file
readr::write_tsv(final_mapped_df, file.path("AnalysisFolderNew",
  results_dir,
  "SRP159810_Symbol_IDs.tsv" # Replace with a relevant output file name
))












