#!/bin/bash
set -euo pipefail

# DNAnexus applet entrypoint for IMPACT-prioritization
# Inputs: $genelist (GeneList.txt), $gds_files (array of .gds files)

# download inputs to working directory
dx download "$genelist"   -o GeneList.txt
for file in "${gds_files[@]}"; do
  dx download "$file" -o .
done

# Confirm R version
R --version

# Run the prioritization script
Rscript IMPACT-prioritization.r

# List files in out/ directory for debugging
if [ -d out ]; then
  echo "Listing files in out/ directory:"
  ls -lh out/
else
  echo "No out/ directory found"
fi

# Prepare output directories
mkdir -p out/snv_impact_gds
mkdir -p out/log_file

# Move files to their respective output directories
mv out/*_SNV_IMPACT.gds out/snv_impact_gds/
mv out/impact_prioritization.log out/log_file/

# Upload all output files and register them
echo "Uploading and registering all output files from out/ directory..."
dx-upload-all-outputs
echo "All output files uploaded and registered."
