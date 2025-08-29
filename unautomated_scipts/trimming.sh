# comment out memory tracking codes if not needed
LOG_FILE="trimming_memory_log.csv"

input_files=(
  "sample1.fastq.gz"
  "sample2.fastq.gz"
  "sample3.fastq.gz"
  "sample4.fastq.gz"
  "sample5.fastq.gz"
  "sample6.fastq.gz"
  "sample7.fastq.gz"
  "sample8.fastq.gz"
  "sample9.fastq.gz"
  "sample10.fastq.gz"
)

for file in "${input_files[@]}"; do
    base=$(basename "$file" .fastq.gz)
    echo "Processing $base..."

    (
        cutadapt -a AGATCGGAAG --cores 8 -o "${base}_trimmed.fastq.gz" "${base}.fastq.gz"
    ) &

    PIPE_PID=$!

    ./memtrack_nopeak.sh $PIPE_PID "$base" "$LOG_FILE" &

    wait $PIPE_PID
done
