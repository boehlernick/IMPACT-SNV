#!/bin/bash
set -euo pipefail

# DNAnexus applet entrypoint for IMPACT-prioritization
# Inputs: $genelist (GeneList.txt), $gds_files (array of .gds files)

# Copy inputs to working directory
dx-download-all-inputs

# Move all downloaded inputs into the working directory
find . -type f -not -path "./*" -exec mv {} . \;

# Confirm R version
R --version

# Run the prioritization script
Rscript IMPACT-prioritization.r

# Upload all output files and register them
echo "Uploading and registering all output files..."
dx-upload-all-outputs
echo "All output files uploaded and registered."
