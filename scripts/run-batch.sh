#!/bin/bash
# scripts/run-batch.sh

TASKS_FILE="/tmp/subsampling_tasks.txt"
MAX_JOBS=3
WAIT_ONLY=0
LOG_DIR="docs/logs/subsampling"

# Generate task list if missing (self-contained)
if [ ! -f "$TASKS_FILE" ]; then
  echo "Generating task list at $TASKS_FILE..."
  : > "$TASKS_FILE"
  for n in 5 10 20 40 70; do
    for r in $(seq 1 50); do
      echo "$n $r" >> "$TASKS_FILE"
    done
  done
  echo "113 1" >> "$TASKS_FILE"
fi

mkdir -p "$LOG_DIR"

trap 'echo -e "\n[!] Timeout reached. Finishing active fits..." ; WAIT_ONLY=1' SIGTERM

if [ ! -f "$TASKS_FILE" ]; then
  echo "Error: $TASKS_FILE not found."
  exit 1
fi

echo "Starting batch run..."

while read -r n r; do
  if [ "$WAIT_ONLY" -eq 1 ]; then break; fi

  while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do wait -n; done
  if [ "$WAIT_ONLY" -eq 1 ]; then break; fi

  out_file=$(printf "results/subsampling/n%03d_rep%02d.rds" "$n" "$r")
  log_file=$(printf "%s/n%03d_rep%02d.log" "$LOG_DIR" "$n" "$r")
  
  if [ -f "$out_file" ]; then
    echo "[SKIP] n=$n rep=$r"
  else
    (
      echo "[RUN]  n=$n rep=$r"
      docker run --rm -e OPENBLAS_NUM_THREADS=1 --memory=4g \
        -v "$PWD":/home/rstudio/project paired-rnaseq:dev \
        Rscript R/09-subsampling.R "$n" "$r" > "$log_file" 2>&1
      
      if [ $? -eq 0 ]; then echo "[OK]   n=$n rep=$r"
      else echo "[FAIL] n=$n rep=$r (check $log_file)"; fi
    ) &
  fi
done < "$TASKS_FILE"

if [ "$WAIT_ONLY" -eq 1 ]; then echo "Waiting for active containers to finish..."; fi
wait
echo "Batch complete."