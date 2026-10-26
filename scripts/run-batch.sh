#!/bin/bash
# scripts/run-batch.sh

TASKS_FILE="/tmp/subsampling_tasks.txt"
MAX_JOBS=3
LOG_DIR="docs/logs/subsampling"

# To do list
if [ ! -f "$TASKS_FILE" ]; then
  echo "Generating task list at $TASKS_FILE..."
  : > "$TASKS_FILE"
  for n in 5 10 20 40 70; do
    for r in $(seq 1 50); do echo "$n $r" >> "$TASKS_FILE"; done
  done
  echo "113 1" >> "$TASKS_FILE"
fi

mkdir -p "$LOG_DIR"
echo "Starting batch run with $MAX_JOBS parallel jobs..."

# Parallel run with xargs
cat "$TASKS_FILE" | xargs -P $MAX_JOBS -n 2 bash -c '
  n=$1
  r=$2
  out_file=$(printf "results/subsampling/n%03d_rep%02d.rds" "$n" "$r")
  
  if [ -f "$out_file" ]; then
    echo "[SKIP] n=$n rep=$r"
  else
    echo "[RUN]  n=$n rep=$r"
    tmp_log=$(mktemp)
    docker run --rm -e OPENBLAS_NUM_THREADS=1 --memory=4g \
      -v "$PWD":/home/rstudio/project paired-rnaseq:dev \
      Rscript R/09-subsampling.R "$n" "$r" > "$tmp_log" 2>&1
    rc=$?
    
    if [ $rc -ne 0 ]; then
      fail_log=$(printf "docs/logs/subsampling/FAIL_n%03d_rep%02d.log" "$n" "$r")
      mv "$tmp_log" "$fail_log"
      echo "[FAIL] n=$n rep=$r (log saved to $fail_log)"
    else
      rm "$tmp_log"
      echo "[OK]   n=$n rep=$r"
    fi
  fi
' _

echo "Batch complete."