
# ============================================================================
# Makefile: Smoke test for IMPACT-SNV Step 1 container
# Requires: Docker installed and Dockerfile + entrypoint.sh present in project root.
# ============================================================================
IMAGE?=impact-step1:latest
SMOKE_DIR:=smoke
DATA_DIR:=$(SMOKE_DIR)/data
OUT_DIR:=$(SMOKE_DIR)/out
PREFIX:=smoke
THREADS?=2

.PHONY: all build-step1 test-step1 clean

all: test-step1

build-step1:
	docker build -t $(IMAGE) .

# Run the container on tiny synthetic VCFs and validate outputs
# Notes:
#  - We mount the entire smoke directory at /work inside the container.
#  - The entrypoint will bgzip/index inputs and produce outputs under /work/out.
#  - We pass an explicit output prefix to make assertions deterministic.
#  - We use --input-vcfs for determinism (predefined list, not directory scan).
test-step1: build-step1
	mkdir -p $(OUT_DIR)
	docker run --rm \
	  -v $(PWD)/$(SMOKE_DIR):/work \
	  $(IMAGE) \
	  --input-vcfs /work/data/NA12892.subset.sorted.vcf.gz,/work/data/NA18591.subset.sorted.vcf.gz,/work/data/NA18596.subset.sorted.vcf.gz \
	  --outdir /work/out \
	  --output-prefix $(PREFIX) \
	  --threads $(THREADS)
	bash $(SMOKE_DIR)/assert.sh

clean:
	rm -rf $(OUT_DIR)
