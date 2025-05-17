### UKB RAP FAVOR Full Database V3.1 ###
### Dependencies ###
# Install Rust and XSV in terminal before running #
# curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# source $HOME/.cargo/env
# rustc --version
# cargo install xsv
# xsv --version

### Usage ###
# Rscript favorannotator_UKB_v2.r FA_UKB_V2_chr15 merged_chr15.gds 15 TRUE
# Note four mandatory inputs in order are "Output filename prefix" "input filename" "chromosome number" "Compression"
# One optional input for the FAVOR CSV database download directly from your own project is the 5th argument

# Start timer
start_time <- Sys.time()

# Get the current PATH
current_path <- Sys.getenv("PATH")

# Append the xsv directory to the PATH
Sys.setenv(PATH = paste(current_path, "/home/dnanexus/.cargo/bin", sep = ":"))



# Optionally, run an xsv command
system("xsv --version")


# Check if xsv is available
xsv_path <- system("which xsv", intern = TRUE)
if (length(xsv_path) == 0) {
  stop("Error: xsv binary not found. Please install xsv or ensure it's accessible in your PATH environment variable.")
}
print(xsv_path)

xsv <- xsv_path

# #Install additional packages
 if (!requireNamespace("BiocManager", quietly=TRUE))
   install.packages("BiocManager")
 BiocManager::install("SeqArray")
 BiocManager::install("SeqVarTools")
 BiocManager::install("gdsfmt")

options(repos = c(CRAN = "http://cran.rstudio.com/"))

install.packages("readr")
install.packages("doParallel")

### R packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(readr)  # For CSV reading
library(parallel)
library(foreach)
library(doParallel)

args <- commandArgs(TRUE)
### mandatory
outfile <- args[1]
gds.file <- args[2]
chr <- as.numeric(args[3])
use_compression <- args[4]

### optional: path to the .txt file with chromosome and File_ID mappings
favor_csv_database <- ifelse(length(args) >= 5, args[5], NULL)

##########################################################################
### Step 0 (Download FAVOR Full Database)
##########################################################################
URLs <- data.frame(chr = c(1:22),
                   URL = c("https://dataverse.harvard.edu/api/access/datafile/6380374",  #1
                           "https://dataverse.harvard.edu/api/access/datafile/6380471",  #2
                           "https://dataverse.harvard.edu/api/access/datafile/6380732",  #3
                           "https://dataverse.harvard.edu/api/access/datafile/6381512",  #4
                           "https://dataverse.harvard.edu/api/access/datafile/6381457",  #5
                           "https://dataverse.harvard.edu/api/access/datafile/6381327",  #6
                           "https://dataverse.harvard.edu/api/access/datafile/6384125",  #7
                           "https://dataverse.harvard.edu/api/access/datafile/6382573",  #8
                           "https://dataverse.harvard.edu/api/access/datafile/6384268",  #9
                           "https://dataverse.harvard.edu/api/access/datafile/6380273",  #10
                           "https://dataverse.harvard.edu/api/access/datafile/6384154",  #11
                           "https://dataverse.harvard.edu/api/access/datafile/6384198",  #12
                           "https://dataverse.harvard.edu/api/access/datafile/6388366",  #13
                           "https://dataverse.harvard.edu/api/access/datafile/6388406",  #14
                           "https://dataverse.harvard.edu/api/access/datafile/6388427",  #15
                           "https://dataverse.harvard.edu/api/access/datafile/6388551",  #16
                           "https://dataverse.harvard.edu/api/access/datafile/6388894",  #17
                           "https://dataverse.harvard.edu/api/access/datafile/6376523",  #18
                           "https://dataverse.harvard.edu/api/access/datafile/6376522",  #19
                           "https://dataverse.harvard.edu/api/access/datafile/6376521",  #20
                           "https://dataverse.harvard.edu/api/access/datafile/6358305",  #21
                           "https://dataverse.harvard.edu/api/access/datafile/6358299")) #22

URL <- URLs[chr, "URL"]
file_name <- gsub(".*?([0-9]+).*", "\\1", URL)  # Extract file name from URL

# Define the first extracted file (e.g., chr21_1.csv for chromosome 21)
first_extracted_file <- paste0("chr", chr, "_1.csv")

# Check if the file already exists in the working directory
if (!file.exists(file_name)) {
  cat("File not found, downloading...\n")
  system(paste0("wget --timeout=60 --progress=bar:force:noscroll ", URL))
} else {
  cat("File already exists, skipping download.\n")
}

# Check if the first extracted file exists to determine if extraction is necessary
if (!file.exists(first_extracted_file)) {
  cat("Extracting file...\n")
  system(paste0("tar -xvf ", file_name))
} else {
  cat("File already extracted, skipping extraction.\n")
}



##########################################################################
### Step 1 (Varinfo_gds)
##########################################################################

### output
output_path <- "./"

### make directory
dir_path <- paste0(output_path, "chr", chr)

# Check if the directory exists
if (dir.exists(dir_path)) {
  # Remove the existing directory and its contents
  unlink(dir_path, recursive = TRUE)
}

# Create the new directory
dir.create(dir_path)


### Chromosome number
## Read info
DB_info <- read.csv(url("https://raw.githubusercontent.com/xihaoli/STAARpipeline-Tutorial/main/FAVORannotator_csv/FAVORdatabase_chrsplit.csv"), header = TRUE)
DB_info <- DB_info[DB_info$Chr == chr, ]

## Open GDS
genofile <- seqOpen(gds.file)

# Print original genofile structure
genofile

CHR <- as.numeric(seqGetData(genofile, "chromosome"))
position <- as.integer(seqGetData(genofile, "position"))
REF <- as.character(seqGetData(genofile, "$ref"))
ALT <- as.character(seqGetData(genofile, "$alt"))

VarInfo_genome <- paste0(CHR, "-", position, "-", REF, "-", ALT)

seqClose(genofile)


## Generate VarInfo
for (kk in 1:nrow(DB_info)) {
  message("Processing chunk: ", kk)
  
  VarInfo <- VarInfo_genome[(position >= DB_info$Start_Pos[kk]) & (position <= DB_info$End_Pos[kk])]
  VarInfo <- data.frame(VarInfo)
  
  write.csv(VarInfo, paste0(output_path, "chr", chr, "/VarInfo_chr", chr, "_", kk, ".csv"), quote = FALSE, row.names = FALSE)
}

##########################################################################
### Step 2 (Annotate)
##########################################################################

### xsv directory
# xsv <- xsv_path

### DB file
DB_path <- "./"

### Anno channel (subset)
anno_colnum <- c(1:190)

chr_splitnum <- sum(DB_info$Chr == chr)


#### Optimization test #####

library(doParallel)

# Set up the parallel backend
cl <- makeCluster(detectCores() - 1) 
registerDoParallel(cl)

# Parallelize annotation
results <- foreach(kk = 1:chr_splitnum, .packages = 'foreach') %dopar% {
  cat("Annotating chunk: ", kk, " for chromosome ", chr, "\n") # Output message to console
  
  # Run system command and redirect output to the expected file
  output_file <- paste0(output_path, "chr", chr, "/Anno_chr", chr, "_", kk, ".csv")
  system_command <- paste0(xsv, " join --left VarInfo ", output_path, "chr", chr, "/VarInfo_chr", chr, "_", kk, ".csv variant_vcf ", 
                           DB_path, "/chr", chr, "_", kk, ".csv > ", output_file, " 2>&1")
  
  # Capture exit status
  exit_status <- system(system_command, intern = FALSE)
  
  # Write user-friendly message based on the exit status
  if (exit_status == 0) {
    paste("Annotating chunk", kk, "for chromosome", chr, "was successful")
  } else {
    paste("Annotating chunk", kk, "for chromosome", chr, "failed with exit code", exit_status)
  }
}

# Print the messages from each chunk
cat(unlist(results), sep = "\n")

# Stop the cluster
stopCluster(cl)




## Merge info
Anno <- paste0(output_path, "chr", chr, "/Anno_chr", chr, "_", seq(1:chr_splitnum), ".csv ")
merge_command <- paste0(xsv, " cat rows ", Anno[1])

for (kk in 2:chr_splitnum) {
  merge_command <- paste0(merge_command, Anno[kk])
}

merge_command <- paste0(merge_command, "> ", output_path, "chr", chr, "/Anno_chr", chr, ".csv")
system(merge_command)

## Subset
anno_colnum_xsv <- paste0(anno_colnum, collapse = ",")
system(paste0(xsv, " select ", anno_colnum_xsv, " ", output_path, "chr", chr, "/Anno_chr", chr, ".csv > ", output_path, "chr", chr, "/Anno_chr", chr, "_STAARpipeline.csv"))

##########################################################################
### Step 3 (gds2agds)
##########################################################################

### Annotation file
dir_anno <- "./"
anno_file_name_1 <- "Anno_chr"
anno_file_name_2 <- "_STAARpipeline.csv"


### Define column types (Some columns can be misinterpreted for some chromosomes depending on sample input)
column_types <- cols(
  aloft_description = col_character(),
  aloft_value = col_double(),
  metasvm_pred = col_character(), 
  clnsig = col_character(),
  clnsigincl = col_character(),
  clndn = col_character(),
  clndnincl = col_character(),
  clnrevstat = col_character(),
  clndisdb = col_character(),
  clndisdbincl = col_character(),
  geneinfo = col_character(),
  refseq_category = col_character(),
  origin = col_character(),
  refseq_info = col_character(),
  sift_cat = col_character(),          # Column 98
  sift_val = col_double(),             # Column 99
  polyphen_cat = col_character(),      # Column 100
  polyphen_val = col_double(),         # Column 101
  grantham = col_double(),             # Column 151
  polyphen2_hdiv_score = col_double(), # Column 49
  polyphen2_hvar_score = col_double(), # Column 50
  mutation_assessor_score = col_double(), # Column 52
  mutation_taster_score = col_double() # Column 51
)

### Read annotation data with specific column types
message("Reading annotation data for chr", chr)
FunctionalAnnotation <- tryCatch(
  {
    read_csv(paste0(dir_anno, "chr", chr, "/", anno_file_name_1, chr, anno_file_name_2), col_types = column_types)
  },
  error = function(e) {
    message("Error reading annotation data: ", e)
    stop("Failed to read annotation data")
  }
)

# Check for parsing issues
parsing_problems <- problems(FunctionalAnnotation)
print(parsing_problems)

if (nrow(parsing_problems) > 0) {
  message("Parsing problems detected. Logging the issues:")
  print(parsing_problems)  # Log parsing problems for debugging
} else {
  message("No parsing issues detected.")
}

# Check for warnings
warning_check <- warnings()  # Get all warnings from the current session
if (length(warning_check) > 0) {
  message("Warnings detected. Logging:")
  print(warning_check)  # Log parsing problems for debugging
} else {
  message("No warnings detected.")
}

# Convert aloft_description to character (string) if needed again
FunctionalAnnotation$aloft_description <- as.character(FunctionalAnnotation$aloft_description)


## Open GDS
message("Opening GDS file for writing annotation.")
genofile <- seqOpen(gds.file, readonly = FALSE)

Anno.folder <- index.gdsn(genofile, "annotation/info")

# Add FunctionalAnnotation to GDS
if (use_compression == "YES") {
  message("Writing annotations with compression.")
  add.gdsn(Anno.folder, "FunctionalAnnotation", val = FunctionalAnnotation, compress = "LZMA_ra", closezip = TRUE)
} else {
  message("Writing annotations without compression.")
  add.gdsn(Anno.folder, "FunctionalAnnotation", val = FunctionalAnnotation)
}
# Print Genofile Structure
genofile

# Close the GDS genofile
seqClose(genofile)

# Move output GDS file
system(paste0("mv ", gds.file, " ", outfile, ".gds"))

message("Process completed successfully for chromosome ", chr)

# Check for warnings again
warning_check2 <- warnings()  # Get all warnings from the current session
if (length(warning_check) > 0) {
  message("Warnings detected. Logging:")
  print(warning_check2)  # Log parsing problems for debugging
} else {
  message("No warnings detected.")
}

# End timer
end_time <- Sys.time()

# Calculate elapsed time
elapsed_time <- end_time - start_time

# Print the elapsed time
cat("Elapsed time:", elapsed_time, "\n")

# Save to a log file
log_file <- "execution_log.txt"
write(paste("Execution started at:", start_time, "\n",
            "Execution ended at:", end_time, "\n",
            "Elapsed time:", elapsed_time), file = log_file, append = TRUE)

# Print arguments for logging
print(args)

