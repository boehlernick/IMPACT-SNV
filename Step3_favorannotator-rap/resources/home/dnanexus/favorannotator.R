### FAVOR Full Database V3.1 Annotator ###
### Designed for DNAnexus Research Analysis Platform ###
### Dependencies ###
# Install Rust and XSV in terminal before running #
# curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# source $HOME/.cargo/env
# rustc --version
# cargo install xsv
# xsv --version

### Usage ###
# Rscript favorannotator.r FA_chr15 merged_chr15.gds 15 TRUE
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


### Hard-coded column types for FAVOR annotation columns
column_types <- cols(
  VarInfo = col_character(),
  vid = col_double(),
  variant_vcf = col_character(),
  variant_annovar = col_character(),
  chromosome = col_double(),
  start_position = col_double(),
  end_position = col_double(),
  ref_annovar = col_character(),
  alt_annovar = col_character(),
  position = col_character(),
  ref_vcf = col_character(),
  alt_vcf = col_character(),
  aloft_value = col_double(),
  aloft_description = col_character(),
  apc_conservation = col_double(),
  apc_conservation_v2 = col_double(),
  apc_epigenetics_active = col_double(),
  apc_epigenetics = col_double(),
  apc_epigenetics_repressed = col_double(),
  apc_epigenetics_transcription = col_double(),
  apc_local_nucleotide_diversity = col_double(),
  apc_local_nucleotide_diversity_v2 = col_double(),
  apc_local_nucleotide_diversity_v3 = col_double(),
  apc_mappability = col_double(),
  apc_micro_rna = col_double(),
  apc_mutation_density = col_double(),
  apc_protein_function = col_double(),
  apc_protein_function_v2 = col_double(),
  apc_protein_function_v3 = col_double(),
  apc_proximity_to_coding = col_double(),
  apc_proximity_to_coding_v2 = col_double(),
  apc_proximity_to_tsstes = col_double(),
  apc_transcription_factor = col_double(),
  bravo_an = col_character(),
  bravo_af = col_character(),
  filter_status = col_character(),
  cage_enhancer = col_character(),
  cage_promoter = col_character(),
  cage_tc = col_character(),
  clnsig = col_character(),
  clnsigincl = col_character(),
  clndn = col_character(),
  clndnincl = col_character(),
  clnrevstat = col_character(),
  origin = col_character(),
  clndisdb = col_character(),
  clndisdbincl = col_character(),
  geneinfo = col_character(),
  polyphen2_hdiv_score = col_double(),
  polyphen2_hvar_score = col_double(),
  mutation_taster_score = col_double(),
  mutation_assessor_score = col_double(),
  metasvm_pred = col_character(),
  fathmm_xf = col_double(),
  funseq_value = col_double(),
  funseq_description = col_character(),
  genecode_comprehensive_category = col_character(),
  genecode_comprehensive_info = col_character(),
  genecode_comprehensive_exonic_category = col_character(),
  genecode_comprehensive_exonic_info = col_character(),
  genehancer = col_character(),
  af_total = col_double(),
  af_asj_female = col_double(),
  af_eas_female = col_double(),
  af_afr_male = col_double(),
  af_female = col_double(),
  af_fin_male = col_double(),
  af_oth_female = col_double(),
  af_ami = col_double(),
  af_oth = col_double(),
  af_male = col_double(),
  af_ami_female = col_double(),
  af_afr = col_double(),
  af_eas_male = col_double(),
  af_sas = col_double(),
  af_nfe_female = col_double(),
  af_asj_male = col_double(),
  af_raw = col_double(),
  af_oth_male = col_double(),
  af_nfe_male = col_double(),
  af_asj = col_double(),
  af_amr_male = col_double(),
  af_amr_female = col_double(),
  af_sas_female = col_double(),
  af_fin = col_double(),
  af_afr_female = col_double(),
  af_sas_male = col_double(),
  af_amr = col_double(),
  af_nfe = col_double(),
  af_eas = col_double(),
  af_ami_male = col_double(),
  af_fin_female = col_double(),
  linsight = col_double(),
  gc = col_double(),
  cpg = col_double(),
  min_dist_tss = col_double(),
  min_dist_tse = col_double(),
  sift_cat = col_character(),
  sift_val = col_double(),
  polyphen_cat = col_character(),
  polyphen_val = col_double(),
  priphcons = col_double(),
  mamphcons = col_double(),
  verphcons = col_double(),
  priphylop = col_double(),
  mamphylop = col_double(),
  verphylop = col_double(),
  bstatistic = col_double(),
  chmm_e1 = col_double(),
  chmm_e2 = col_double(),
  chmm_e3 = col_double(),
  chmm_e4 = col_double(),
  chmm_e5 = col_double(),
  chmm_e6 = col_double(),
  chmm_e7 = col_double(),
  chmm_e8 = col_double(),
  chmm_e9 = col_double(),
  chmm_e10 = col_double(),
  chmm_e11 = col_double(),
  chmm_e12 = col_double(),
  chmm_e13 = col_double(),
  chmm_e14 = col_double(),
  chmm_e15 = col_double(),
  chmm_e16 = col_double(),
  chmm_e17 = col_double(),
  chmm_e18 = col_double(),
  chmm_e19 = col_double(),
  chmm_e20 = col_double(),
  chmm_e21 = col_double(),
  chmm_e22 = col_double(),
  chmm_e23 = col_double(),
  chmm_e24 = col_double(),
  chmm_e25 = col_double(),
  gerp_rs = col_double(),
  gerp_rs_pval = col_double(),
  gerp_n = col_double(),
  gerp_s = col_double(),
  encodeh3k4me1_sum = col_double(),
  encodeh3k4me2_sum = col_double(),
  encodeh3k4me3_sum = col_double(),
  encodeh3k9ac_sum = col_double(),
  encodeh3k9me3_sum = col_double(),
  encodeh3k27ac_sum = col_double(),
  encodeh3k27me3_sum = col_double(),
  encodeh3k36me3_sum = col_double(),
  encodeh3k79me2_sum = col_double(),
  encodeh4k20me1_sum = col_double(),
  encodeh2afz_sum = col_double(),
  encode_dnase_sum = col_double(),
  encodetotal_rna_sum = col_double(),
  grantham = col_double(),
  freq100bp = col_double(),
  rare100bp = col_double(),
  sngl100bp = col_double(),
  freq1000bp = col_double(),
  rare1000bp = col_double(),
  sngl1000bp = col_double(),
  freq10000bp = col_double(),
  rare10000bp = col_double(),
  sngl10000bp = col_double(),
  remap_overlap_tf = col_double(),
  remap_overlap_cl = col_double(),
  cadd_rawscore = col_double(),
  cadd_phred = col_double(),
  k24_bismap = col_double(),
  k24_umap = col_double(),
  k36_bismap = col_double(),
  k36_umap = col_double(),
  k50_bismap = col_double(),
  k50_umap = col_double(),
  k100_bismap = col_double(),
  k100_umap = col_double(),
  nucdiv = col_double(),
  rdhs = col_character(),
  recombination_rate = col_double(),
  refseq_category = col_character(),
  refseq_info = col_character(),
  refseq_exonic_category = col_character(),
  refseq_exonic_info = col_character(),
  super_enhancer = col_character(),
  tg_afr = col_double(),
  tg_all = col_double(),
  tg_amr = col_double(),
  tg_eas = col_double(),
  tg_eur = col_double(),
  tg_sas = col_double(),
  ucsc_category = col_character(),
  ucsc_info = col_character(),
  ucsc_exonic_category = col_character(),
  ucsc_exonic_info = col_character()
)

# Function to update column_types for problematic columns
update_column_types_to_character <- function(column_types, problem_cols) {
  for (col in problem_cols) {
    column_types[[col]] <- col_character()
  }
  return(column_types)
}

### Read annotation data with specific column types
message("Reading annotation data for chr", chr)
anno_csv_path <- paste0(dir_anno, "chr", chr, "/", anno_file_name_1, chr, anno_file_name_2)

# Initial read
FunctionalAnnotation <- tryCatch(
  {
    read_csv(anno_csv_path, col_types = column_types)
  },
  error = function(e) {
    message("Error reading annotation data: ", e)
    stop("Failed to read annotation data")
  }
)

# Check for parsing issues and fix them automatically
parsing_problems <- problems(FunctionalAnnotation)
if (nrow(parsing_problems) > 0) {
  message("Parsing problems detected. Attempting to fix by converting problematic columns to string:")
  problem_cols <- unique(parsing_problems$col)
  for (col in problem_cols) {
    colname <- names(FunctionalAnnotation)[as.integer(col)]
    warning(sprintf("Column '%s' (index %s) has parsing issues and will be forced to string.", colname, col))
  }
  # Update column_types for problematic columns
  new_column_types <- update_column_types_to_character(column_types, as.integer(problem_cols))
  # Re-read with updated column_types
  FunctionalAnnotation <- read_csv(anno_csv_path, col_types = new_column_types)
  # Check again for parsing problems
  parsing_problems <- problems(FunctionalAnnotation)
  if (nrow(parsing_problems) > 0) {
    stop("Parsing problems remain after forcing problematic columns to string. Please check the input file.")
  }
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


## Open GDS
message("Opening GDS file for writing annotation.")
genofile <- seqOpen(gds.file, readonly = FALSE)

Anno.folder <- index.gdsn(genofile, "annotation/info")

# Before adding FunctionalAnnotation to GDS, check if node exists and delete if so
if ("FunctionalAnnotation" %in% ls.gdsn(Anno.folder)) {
  delete.gdsn(index.gdsn(Anno.folder, "FunctionalAnnotation"))
}

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

