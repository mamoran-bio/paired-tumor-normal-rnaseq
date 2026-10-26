#!/bin/bash
# scripts/run-loo.sh

MAX_JOBS=3
LOG_DIR="docs/logs/loo"
mkdir -p "$LOG_DIR"

# Get the list of 113 patients (from the model-inputs.rds object)
# We use a small R command to extract them
echo "Extracting patient list..."
docker run --rm -v "$PWD":/home/rstudio/project paired-rnaseq:dev \
  Rscript -e 'inputs <- readRDS("results/model-inputs.rds"); cat(unique(inputs$coldata$patient), sep="\n")' > /tmp/loo_tasks.txt

TOTAL=$(wc -l < /tmp/loo_tasks.txt)
echo "Starting Leave-One-Out for $TOTAL patients with $MAX_JOBS parallel jobs..."

cat /tmp/loo_tasks.txt | xargs -P $MAX_JOBS -I {} bash -c '
  patient=$1
  out_file="results/loo/loo_${patient}.rds"
  
  if [ -f "$out_file" ]; then
    echo "[SKIP] $patient"
  else
    echo "[RUN]  $patient"
    tmp_log=$(mktemp)
    docker run --rm -e OPENBLAS_NUM_THREADS=1 --memory=4g \
      -v "$PWD":/home/rstudio/project paired-rnaseq:dev \
      Rscript R/10-leave-one-out.R "$patient" > "$tmp_log" 2>&1
    rc=$?
    
    if [ $rc -ne 0 ]; then
      fail_log=$(printf "docs/logs/loo/FAIL_%s.log" "$patient")
      mv "$tmp_log" "$fail_log"
      echo "[FAIL] $patient (log saved to $fail_log)"
    else
      rm "$tmp_log"
      echo "[OK]   $patient"
    fi
  fi
' _ {}

echo "Leave-One-Out complete."