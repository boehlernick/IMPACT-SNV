#!/usr/bin/env bash
set -euo pipefail

ts() { date +"[%Y-%m-%d %H:%M:%S]"; }
die() { echo "$(ts) ERROR: $*" >&2; exit 1; }
have_cmd() { command -v "$1" >/dev/null 2>&1; }

usage() {
  cat <<'USAGE'
IMPACT-SNV Step 1: VCF merge + normalization + per-contig split
================================================================
Usage:
  docker run --rm -v $PWD:/work impact-step1 \
    --input-vcfs /work/a.vcf.gz,/work/b.vcf.gz \
    --outdir /work/out \
    [--threads 8] \
    [--output-prefix merged_output] \
    [--keep-format] \
    [--keep-format-fields DP,GQ]

  OR

  docker run --rm -v $PWD:/work impact-step1 \
    --input-dir /work/vcfs \
    --outdir /work/out \
    [--threads 8] \
    [--output-prefix merged_output] \
    [--keep-format-fields DP,GQ]

Input Modes (choose one):
  --input-vcfs          Comma-separated list of VCF paths (bgzipped or plain)
  --input-dir           Directory containing .vcf or .vcf.gz files

FORMAT Processing (mutually exclusive):
  (default)             Strip all FORMAT fields except GT
  --keep-format         Preserve all FORMAT fields from inputs
  --keep-format-fields  Keep specific FORMAT fields (comma-separated list)
                        GT is always included; this adds to it

Optional:
  --output-prefix       Base name for merged outputs (default: merged_output)
  --outdir              Output directory (default: /out)
  --threads             Parallel threads for bgzip (default: 4)

Outputs:
  <outdir>/<prefix>.vcf.gz + .tbi
  <outdir>/<prefix>.normalized.vcf.gz + .tbi
  <outdir>/<prefix>.normalized.Chr1.vcf.gz + .tbi
  ... (20 more per-chromosome files)
  <outdir>/<prefix>.normalized.ChrX.vcf.gz + .tbi
  <outdir>/<prefix>.normalized.ChrY.vcf.gz + .tbi
  <outdir>/<prefix>.normalized.ChrMT.vcf.gz + .tbi

Total: 25 contig files + 25 indexes (50 files).
Contigs normalized to canonical names: Chr1, Chr2, ..., Chr22, ChrX, ChrY, ChrMT.
Missing contigs receive empty header-only files.

Chromosome Naming:
  Input aliases (1/chr1/Chr1, etc.) automatically map to canonical names.
  Output contig names are standardized across all 25 files.
USAGE
}

# ============================================================================
# Canonical contigs and their aliases for bcftools view --regions
# ============================================================================
declare -A CANONICAL_CONTIGS=(
  [Chr1]="Chr1,chr1,1"
  [Chr2]="Chr2,chr2,2"
  [Chr3]="Chr3,chr3,3"
  [Chr4]="Chr4,chr4,4"
  [Chr5]="Chr5,chr5,5"
  [Chr6]="Chr6,chr6,6"
  [Chr7]="Chr7,chr7,7"
  [Chr8]="Chr8,chr8,8"
  [Chr9]="Chr9,chr9,9"
  [Chr10]="Chr10,chr10,10"
  [Chr11]="Chr11,chr11,11"
  [Chr12]="Chr12,chr12,12"
  [Chr13]="Chr13,chr13,13"
  [Chr14]="Chr14,chr14,14"
  [Chr15]="Chr15,chr15,15"
  [Chr16]="Chr16,chr16,16"
  [Chr17]="Chr17,chr17,17"
  [Chr18]="Chr18,chr18,18"
  [Chr19]="Chr19,chr19,19"
  [Chr20]="Chr20,chr20,20"
  [Chr21]="Chr21,chr21,21"
  [Chr22]="Chr22,chr22,22"
  [ChrX]="ChrX,chrX,X"
  [ChrY]="ChrY,chrY,Y"
  [ChrMT]="ChrMT,chrMT,MT,chrM,M"
)

# List of canonical contigs in order
CANONICAL_ORDER=(Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 Chr9 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 Chr20 Chr21 Chr22 ChrX ChrY ChrMT)

# ============================================================================
# Helper function: get aliases for a canonical contig
# ============================================================================
aliases_for() {
  local contig="$1"
  echo "${CANONICAL_CONTIGS[$contig]}"
}

# ============================================================================
# Helper function: normalize contig name to canonical form
# ============================================================================
# ============================================================================
# Helper function: normalize contig name to canonical form (case-safe, exact match)
# ============================================================================
normalize_contig_name() {
  local name="$1"
  local name_upper="${name^^}"  # Convert to uppercase
  
  # Try to find which canonical contig this name belongs to
  for canonical in "${!CANONICAL_CONTIGS[@]}"; do
    local aliases="${CANONICAL_CONTIGS[$canonical]}"
    # Split aliases by comma and check each one (case-safe, exact match)
    IFS=',' read -ra alias_tokens <<< "$aliases"
    for token in "${alias_tokens[@]}"; do
      local token_upper="${token^^}"  # Convert token to uppercase
      if [[ "$name_upper" == "$token_upper" ]]; then
        echo "$canonical"
        return 0
      fi
    done
  done
  
  # If not found, return as-is (will fail downstream gracefully)
  echo "$name"
}

# ============================================================================
# Helper function: create empty header-only VCF for missing contigs
# ============================================================================
create_empty_contig_vcf() {
  local header_vcf="$1"
  local contig="$2"
  local output="$3"
  local tmpheader="${output}.tmp.header"
  local tmpvcf="${output}.tmp.vcf.gz"
  
  echo "$(ts) INFO: Creating empty VCF for $contig: $output"
  
  # Extract header from normalized VCF
  bcftools view -h "$header_vcf" > "$tmpheader"
  
  # Append contig line if not already present (use bcftools reheader for proper handling)
  if ! grep -Fq "##contig=<ID=${contig}" "$tmpheader"; then
    # Try to find the contig line in the source file using case-insensitive matching
    # Map canonical names to possible source names (e.g., Chr1 -> chr1, ChrMT -> chrM)
    local source_contig="${contig}"
    case "$contig" in
      Chr[0-9]*|ChrX|ChrY) source_contig=$(echo "$contig" | tr '[:upper:]' '[:lower:]') ;;
      ChrMT) source_contig="chrM" ;;
    esac
    
    # Try to extract the contig line from source, or use a reasonable default
    local contig_line=$(bcftools view -h "$header_vcf" | grep "^##contig=<ID=${source_contig}" || echo "")
    if [[ -z "$contig_line" ]]; then
      # Fallback: use reasonable defaults based on GRCh38
      case "$contig" in
        Chr1) contig_line="##contig=<ID=${contig},length=248956422>" ;;
        Chr2) contig_line="##contig=<ID=${contig},length=242193529>" ;;
        Chr3) contig_line="##contig=<ID=${contig},length=198295559>" ;;
        Chr4) contig_line="##contig=<ID=${contig},length=190214555>" ;;
        Chr5) contig_line="##contig=<ID=${contig},length=181538259>" ;;
        Chr6) contig_line="##contig=<ID=${contig},length=170805979>" ;;
        Chr7) contig_line="##contig=<ID=${contig},length=159345973>" ;;
        Chr8) contig_line="##contig=<ID=${contig},length=145138636>" ;;
        Chr9) contig_line="##contig=<ID=${contig},length=138394717>" ;;
        Chr10) contig_line="##contig=<ID=${contig},length=133797422>" ;;
        Chr11) contig_line="##contig=<ID=${contig},length=135086622>" ;;
        Chr12) contig_line="##contig=<ID=${contig},length=133275309>" ;;
        Chr13) contig_line="##contig=<ID=${contig},length=114364328>" ;;
        Chr14) contig_line="##contig=<ID=${contig},length=107043718>" ;;
        Chr15) contig_line="##contig=<ID=${contig},length=101991189>" ;;
        Chr16) contig_line="##contig=<ID=${contig},length=90338345>" ;;
        Chr17) contig_line="##contig=<ID=${contig},length=83257441>" ;;
        Chr18) contig_line="##contig=<ID=${contig},length=80373285>" ;;
        Chr19) contig_line="##contig=<ID=${contig},length=58617616>" ;;
        Chr20) contig_line="##contig=<ID=${contig},length=64444167>" ;;
        Chr21) contig_line="##contig=<ID=${contig},length=46709983>" ;;
        Chr22) contig_line="##contig=<ID=${contig},length=50818468>" ;;
        ChrX) contig_line="##contig=<ID=${contig},length=156040895>" ;;
        ChrY) contig_line="##contig=<ID=${contig},length=57227415>" ;;
        ChrMT) contig_line="##contig=<ID=${contig},length=16569>" ;;
        *) contig_line="##contig=<ID=${contig}>" ;;
      esac
    else
      # Use the found line, but replace the source ID with our canonical ID
      contig_line=$(echo "$contig_line" | sed "s/ID=[^,]*/ID=${contig}/")
    fi
    echo "$contig_line" >> "$tmpheader"
  fi
  
  # Create empty VCF by filtering to nonexistent region, then reheader it
  bcftools view "$header_vcf" -r "nonexistent_region_xyz_12345" -Oz -o "$tmpvcf" 2>/dev/null || true
  
  # Use bcftools reheader to apply the modified header (this preserves our added contig lines)
  bcftools reheader -h "$tmpheader" "$tmpvcf" -o "$output"
  
  # Index the output
  bcftools index -f "$output"
  
  # Clean up temp files
  rm -f "$tmpheader" "$tmpvcf" "${tmpvcf}.csi" "${tmpvcf}.tbi"
}

# ============================================================================
# Parse arguments and set defaults
# ============================================================================
INPUT_VCFS=""
INPUT_DIR=""
OUTDIR="/out"
OUTPUT_PREFIX="merged_output"
THREADS="${THREADS:-4}"
KEEP_FORMAT="false"
KEEP_FORMAT_FIELDS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-vcfs) INPUT_VCFS="$2"; shift 2;;
    --input-dir) INPUT_DIR="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --output-prefix) OUTPUT_PREFIX="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --keep-format) KEEP_FORMAT="true"; shift 1;;
    --keep-format-fields) KEEP_FORMAT_FIELDS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) die "Unknown argument: $1";;
  esac
done

# Validate input mode
if [[ -z "$INPUT_VCFS" && -z "$INPUT_DIR" ]]; then
  usage
  die "Either --input-vcfs or --input-dir is required"
fi

if [[ -n "$INPUT_VCFS" && -n "$INPUT_DIR" ]]; then
  die "Cannot use both --input-vcfs and --input-dir; choose one"
fi

if [[ "$KEEP_FORMAT" == "true" && -n "$KEEP_FORMAT_FIELDS" ]]; then
  die "Cannot use both --keep-format and --keep-format-fields"
fi

mkdir -p "$OUTDIR"

# Verify required tools
for cmd in bcftools tabix bgzip; do
  have_cmd "$cmd" || die "Missing required tool: $cmd"
done

# ============================================================================
# Phase 1: Input Discovery & Preparation
# ============================================================================

echo "$(ts) INFO: Discovering and preparing input VCFs"

CLEANED_VCFS=()

if [[ -n "$INPUT_VCFS" ]]; then
  # Parse comma-separated list
  IFS=',' read -r -a VCF_ARR <<< "$INPUT_VCFS"
  for v in "${VCF_ARR[@]}"; do
    [[ -f "$v" ]] || die "Input file not found: $v"
  done
  INPUT_FILES=("${VCF_ARR[@]}")
elif [[ -n "$INPUT_DIR" ]]; then
  # Discover .vcf and .vcf.gz files
  [[ -d "$INPUT_DIR" ]] || die "Input directory not found: $INPUT_DIR"
  mapfile -t INPUT_FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) | sort)
  [[ ${#INPUT_FILES[@]} -gt 0 ]] || die "No .vcf or .vcf.gz files found in $INPUT_DIR"
  echo "$(ts) INFO: Found ${#INPUT_FILES[@]} VCF file(s) in $INPUT_DIR"
fi

# Process each input file
for v in "${INPUT_FILES[@]}"; do
  echo "$(ts) INFO: Processing $v"
  
  # Ensure it's bgzipped
  gzvcf="$v"
  if [[ "$v" != *.vcf.gz ]]; then
    echo "$(ts) INFO: bgzipping $v"
    bgzip -@ "$THREADS" -f "$v"
    gzvcf="${v}.gz"
  fi
  
  # Index if missing
  if [[ ! -f "${gzvcf}.tbi" && ! -f "${gzvcf}.csi" ]]; then
    echo "$(ts) INFO: Indexing $gzvcf"
    bcftools index -f "$gzvcf"
  fi
  
  # Apply FORMAT filtering
  if [[ "$KEEP_FORMAT" == "false" && -z "$KEEP_FORMAT_FIELDS" ]]; then
    # Default: keep only GT
    echo "$(ts) INFO: Stripping FORMAT fields (keeping GT only) for $v"
    tmp="${OUTDIR}/.tmp.$(basename "${gzvcf%%.vcf.gz}").gtonly.vcf.gz"
    bcftools annotate -x 'FORMAT,^GT' "$gzvcf" -Oz -o "$tmp"
    bcftools index -f "$tmp"
    CLEANED_VCFS+=("$tmp")
  elif [[ -n "$KEEP_FORMAT_FIELDS" ]]; then
    # Keep specified fields plus GT
    echo "$(ts) INFO: Keeping FORMAT fields: GT,$KEEP_FORMAT_FIELDS for $v"
    tmp="${OUTDIR}/.tmp.$(basename "${gzvcf%%.vcf.gz}").filtered.vcf.gz"
    # Build the FORMAT expression: GT and the user-specified fields
    format_expr="FORMAT,^GT"
    IFS=',' read -r -a field_arr <<< "$KEEP_FORMAT_FIELDS"
    for field in "${field_arr[@]}"; do
      format_expr="$format_expr,$field"
    done
    bcftools annotate -x "$format_expr" "$gzvcf" -Oz -o "$tmp"
    bcftools index -f "$tmp"
    CLEANED_VCFS+=("$tmp")
  else
    # --keep-format: preserve all
    echo "$(ts) INFO: Preserving all FORMAT fields for $v"
    CLEANED_VCFS+=("$gzvcf")
  fi
done

[[ ${#CLEANED_VCFS[@]} -gt 0 ]] || die "No valid VCF files to merge"

# ============================================================================
# Phase 2: Merge
# ============================================================================

MERGED="${OUTDIR}/${OUTPUT_PREFIX}.vcf.gz"
echo "$(ts) INFO: Merging ${#CLEANED_VCFS[@]} VCF file(s) -> $MERGED"
bcftools merge -O z -o "$MERGED" "${CLEANED_VCFS[@]}"
bcftools index -f "$MERGED"

# ============================================================================
# Phase 3: Normalize & Output
# ============================================================================

NORM="${OUTDIR}/${OUTPUT_PREFIX}.normalized.vcf.gz"
echo "$(ts) INFO: Normalizing multiallelic variants -> $NORM"
bcftools norm -m -any "$MERGED" -Oz -o "$NORM"
bcftools index -f "$NORM"

# ============================================================================
# Phase 4: Split by Canonical Contig & Create Empty VCFs
# ============================================================================

echo "$(ts) INFO: Splitting into 25 canonical contigs"

# Track which canonical contigs exist in the normalized VCF
declare -A contigs_found

# Get list of unique contigs in normalized VCF
mapfile -t actual_contigs < <(bcftools query -f '%CHROM\n' "$NORM" | sort -u)
for ac in "${actual_contigs[@]}"; do
  canonical=$(normalize_contig_name "$ac")
  contigs_found[$canonical]=1
done

# Create per-contig VCFs for all canonical contigs
for contig in "${CANONICAL_ORDER[@]}"; do
  out="${OUTDIR}/${OUTPUT_PREFIX}.normalized.${contig}.vcf.gz"
  
  if [[ -n ${contigs_found[$contig]+x} ]]; then
    # Contig exists; extract it with aliases
    aliases=$(aliases_for "$contig")
    echo "$(ts) INFO: Extracting $contig (aliases: $aliases) -> $out"
    bcftools view "$NORM" --regions "$aliases" -Oz -o "$out"
    bcftools index -f "$out"
  else
    # Contig missing; create empty header-only VCF
    create_empty_contig_vcf "$NORM" "$contig" "$out"
  fi
done

# ============================================================================
# Phase 5: Cleanup & Summary
# ============================================================================

# Remove temporary files
echo "$(ts) INFO: Cleaning up temporary files"
rm -f "$OUTDIR"/.tmp.*.vcf.gz "$OUTDIR"/.tmp.*.vcf.gz.tbi "$OUTDIR"/.tmp.*.vcf.gz.csi

# List outputs
echo "$(ts) INFO: Done. Outputs in $OUTDIR:"
echo "$(ts) INFO: - ${OUTPUT_PREFIX}.vcf.gz + .tbi/.csi"
echo "$(ts) INFO: - ${OUTPUT_PREFIX}.normalized.vcf.gz + .tbi/.csi"
echo "$(ts) INFO: - 25 per-contig VCFs (${OUTPUT_PREFIX}.normalized.Chr1.vcf.gz ... ChrMT.vcf.gz) + indexes"
echo "$(ts) INFO: Total: 54 files (27 VCFs + 27 indexes)"

# Final sanity check
expected_vcfs=27  # merged + normalized + 25 contigs
expected_indexes=27  # merged + normalized + 25 contigs (all get indexes)
actual_vcfs=$(find "$OUTDIR" -name "${OUTPUT_PREFIX}*.vcf.gz" -type f | wc -l)
actual_indexes=$(find "$OUTDIR" \( -name "${OUTPUT_PREFIX}*.vcf.gz.tbi" -o -name "${OUTPUT_PREFIX}*.vcf.gz.csi" \) -type f | wc -l)

echo "$(ts) INFO: VCF files: $actual_vcfs (expected $expected_vcfs)"
echo "$(ts) INFO: Index files: $actual_indexes (expected $expected_indexes)"

if [[ $actual_vcfs -lt $expected_vcfs ]]; then
  die "Expected at least $expected_vcfs VCF files, found $actual_vcfs"
fi

if [[ $actual_indexes -lt $expected_indexes ]]; then
  die "Expected at least $expected_indexes index files, found $actual_indexes"
fi

echo "$(ts) INFO: Step 1 complete!"
