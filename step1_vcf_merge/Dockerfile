
FROM ubuntu:22.04

LABEL org.opencontainers.image.title="IMPACT-SNV Step 1: VCF Merge & Normalization" \
      org.opencontainers.image.description="Merges multi-sample VCFs, normalizes alleles, and splits into 25 canonical contigs" \
      org.opencontainers.image.vendor="IMPACT pipeline" \
      org.opencontainers.image.licenses="GPL-3.0" \
      org.opencontainers.image.version="1.0" \
      maintainer="IMPACT dev team"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      bcftools \
      tabix \
      samtools \
      bash \
      coreutils \
      ca-certificates \
      procps \
      pigz \
      grep \
      findutils && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app
RUN mkdir -p /work /out

COPY entrypoint.sh /app/entrypoint.sh
RUN chmod +x /app/entrypoint.sh

ENV THREADS=4

ENTRYPOINT ["/app/entrypoint.sh"]
CMD ["--help"]
