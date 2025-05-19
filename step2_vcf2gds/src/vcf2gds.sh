#!/bin/bash
# vcf2gds 0.0.2
# Andrew R Wood
# University of Exeter
# No warranty provided!
# App will run on a single machine from beginning to end.
# Modified by Nicholas Boehler, University of Toronto, for the IMPACT pipeline to process an array of VCF files.

main() {
    # Set a default value for $parallel if not specified by the user
    if [ -z "$parallel" ]; then
        parallel=4
        echo "No parallel value specified. Using default: $parallel"
    else
        echo "Number of parallel processes for file conversion: '$parallel'"
    fi

    # Download all VCF files
    for vcf_file in "${vcf_files[@]}"; do
        echo "Downloading VCF: '$vcf_file'"
        dx download "$vcf_file"
    done

    # Create an output directory for GDS files
    mkdir -p /home/dnanexus/out/gds_files

    # Loop through all VCF files in the working directory
    for vcf_path in *.vcf *.vcf.gz; do
        if [ -f "$vcf_path" ]; then
            echo "Processing VCF: '$vcf_path'"

            # Extract the base name for the output GDS file
            base_name=$(basename "$vcf_path" .vcf.gz)
            base_name=$(basename "$base_name" .vcf)
            gds_filename="${base_name}.gds"
            echo "Output GDS: '$gds_filename'"

            # Unpack the R library
            echo "Unpacking R library..."
            tar -zxf r_library.tar.gz
            echo "R library unpacked."

            # Run the R script to 1) convert VCF to GDS, 2) inject PASS as a QC filter
            echo "Running R script with VCF file path: '$vcf_path'..."
            Rscript vcf2gds.R "$vcf_path" "$gds_filename" $parallel
            echo "R script completed."

            # Check if the GDS file was created
            if [ -f "$gds_filename" ]; then
                echo "GDS file created: '$gds_filename'"
                # Move the GDS file to the output directory
                mv "$gds_filename" /home/dnanexus/out/gds_files/
            else
                echo "Error: GDS file not created."
                exit 1
            fi
        fi
    done

    # Upload all output files and register them
    echo "Uploading and registering all output files..."
    dx-upload-all-outputs
    echo "All output files uploaded and registered."
}