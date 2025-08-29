# Memory track - uncomment all the lines to enable memory tracking
#LOG_FILE="kraken_memory_log.csv"

for file in sample*_trimmed.fastq.gz; do
    base=$(basename "$file" _trimmed.fastq.gz)
    echo "Processing $base..."

    kraken2 --db ../kraken2_standard --gzip-compressed --threads 8 --output kraken_${base}.out --use-names --report-minimizer-data  --report kraken_${base}_report "$file" &

    #PIPE_PID=$!

    
    #./memtrack_nopeak.sh $PIPE_PID "$base" "$LOG_FILE" &

    #wait $PIPE_PID
done