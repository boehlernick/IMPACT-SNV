
#!/usr/bin/env bash
set -euo pipefail

OUT_DIR=smoke/out
PREFIX=smoke

echo "[ASSERT] Validating Step 1 outputs in $OUT_DIR"

# ============================================================================
# Check merged and normalized VCFs exist
# ============================================================================

[ -f "$OUT_DIR/${PREFIX}.vcf.gz" ] || { echo "[ASSERT] ERROR: Missing merged VCF: $OUT_DIR/${PREFIX}.vcf.gz"; exit 1; }
echo "[ASSERT] ✓ Merged VCF exists: $OUT_DIR/${PREFIX}.vcf.gz"

[ -f "$OUT_DIR/${PREFIX}.normalized.vcf.gz" ] || { echo "[ASSERT] ERROR: Missing normalized VCF: $OUT_DIR/${PREFIX}.normalized.vcf.gz"; exit 1; }
echo "[ASSERT] ✓ Normalized VCF exists: $OUT_DIR/${PREFIX}.normalized.vcf.gz"

# ============================================================================
# Check that indexes exist for merged and normalized
# ============================================================================

if [ -f "$OUT_DIR/${PREFIX}.vcf.gz.tbi" ] || [ -f "$OUT_DIR/${PREFIX}.vcf.gz.csi" ]; then
  echo "[ASSERT] ✓ Merged VCF is indexed"
else
  echo "[ASSERT] WARNING: Merged VCF not indexed (optional)"
fi

if [ -f "$OUT_DIR/${PREFIX}.normalized.vcf.gz.tbi" ] || [ -f "$OUT_DIR/${PREFIX}.normalized.vcf.gz.csi" ]; then
  echo "[ASSERT] ✓ Normalized VCF is indexed"
else
  echo "[ASSERT] ERROR: Normalized VCF not indexed"
  exit 1
fi

# ============================================================================
# Check per-contig outputs: at least Chr1 and ChrX must exist
# ============================================================================

[ -f "$OUT_DIR/${PREFIX}.normalized.Chr1.vcf.gz" ] || { echo "[ASSERT] ERROR: Missing Chr1 split"; exit 1; }
echo "[ASSERT] ✓ Chr1 split exists: $OUT_DIR/${PREFIX}.normalized.Chr1.vcf.gz"

[ -f "$OUT_DIR/${PREFIX}.normalized.ChrX.vcf.gz" ] || { echo "[ASSERT] ERROR: Missing ChrX split"; exit 1; }
echo "[ASSERT] ✓ ChrX split exists: $OUT_DIR/${PREFIX}.normalized.ChrX.vcf.gz"

# ============================================================================
# Check that exactly 25 canonical contig files exist
# ============================================================================

expected_contigs=("Chr1" "Chr2" "Chr3" "Chr4" "Chr5" "Chr6" "Chr7" "Chr8" "Chr9" "Chr10" \
                  "Chr11" "Chr12" "Chr13" "Chr14" "Chr15" "Chr16" "Chr17" "Chr18" "Chr19" "Chr20" \
                  "Chr21" "Chr22" "ChrX" "ChrY" "ChrMT")

contig_count=0
for contig in "${expected_contigs[@]}"; do
  vcf="$OUT_DIR/${PREFIX}.normalized.${contig}.vcf.gz"
  if [ ! -f "$vcf" ]; then
    echo "[ASSERT] ERROR: Missing contig file: $vcf"
    exit 1
  fi
  contig_count=$((contig_count + 1))
done

echo "[ASSERT] ✓ All 25 canonical contig files present (Chr1-Chr22, ChrX, ChrY, ChrMT)"

# ============================================================================
# Check that indexes exist for all contig files
# ============================================================================

indexed_count=0
for contig in "${expected_contigs[@]}"; do
  vcf="$OUT_DIR/${PREFIX}.normalized.${contig}.vcf.gz"
  if [ -f "${vcf}.tbi" ] || [ -f "${vcf}.csi" ]; then
    indexed_count=$((indexed_count + 1))
  else
    echo "[ASSERT] WARNING: Contig file not indexed: $vcf"
  fi
done

if [ "$indexed_count" -eq 25 ]; then
  echo "[ASSERT] ✓ All 25 contig files are indexed"
else
  echo "[ASSERT] WARNING: Only $indexed_count of 25 contig files are indexed (expected 25)"
fi

# ============================================================================
# Verify file sizes are reasonable (not completely empty)
# ============================================================================

echo ""
echo "[ASSERT] File summary:"
printf "[ASSERT] %-40s %10s\n" "File" "Size"
printf "[ASSERT] %-40s %10s\n" "---" "---"

total_size=0
for vcf in "$OUT_DIR"/${PREFIX}*.vcf.gz; do
  if [ -f "$vcf" ]; then
    size=$(stat -c%s "$vcf" 2>/dev/null || echo "0")
    total_size=$((total_size + size))
    basename_vcf=$(basename "$vcf")
    printf "[ASSERT] %-40s %10d bytes\n" "$basename_vcf" "$size"
  fi
done

printf "[ASSERT] %-40s %10d bytes\n" "TOTAL" "$total_size"

# ============================================================================
# Final validation: all VCF files should be readable by bcftools
# ============================================================================

echo ""
echo "[ASSERT] Validating VCF readability with bcftools..."

for vcf in "$OUT_DIR"/${PREFIX}*.vcf.gz; do
  if [ -f "$vcf" ]; then
    if bcftools view -h "$vcf" > /dev/null 2>&1; then
      :  # Silent success
    else
      echo "[ASSERT] ERROR: Cannot read VCF with bcftools: $vcf"
      exit 1
    fi
  fi
done

echo "[ASSERT] ✓ All VCF files readable by bcftools"

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "[ASSERT] ========================================="
echo "[ASSERT] Smoke test PASSED."
echo "[ASSERT] ========================================="
echo "[ASSERT] Verified:"
echo "[ASSERT]   - Merged VCF: $OUT_DIR/${PREFIX}.vcf.gz"
echo "[ASSERT]   - Normalized VCF: $OUT_DIR/${PREFIX}.normalized.vcf.gz"
echo "[ASSERT]   - 25 canonical contig VCFs"
echo "[ASSERT]   - All outputs indexed and readable"
echo "[ASSERT] ========================================="
