#!/usr/bin/env bash
set -euo pipefail

echo "Checking Python..."
python -c "import sys; print(sys.version)"

echo "Checking R packages..."
Rscript - <<'EOF'
pkgs <- c("optparse", "arrow", "SeqArray", "gdsfmt")
ok <- sapply(pkgs, requireNamespace, quietly = TRUE)
print(ok)
if (!all(ok)) stop("Missing packages")
EOF

echo "✅ Environment OK"
