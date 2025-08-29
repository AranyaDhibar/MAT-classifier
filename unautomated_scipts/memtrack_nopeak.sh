#!/bin/bash

# Usage: ./memtrack_nopeak.sh <PID> <SAMPLE_NAME> <OUTPUT_FILE>
PID=$1
SAMPLE=$2
OUTFILE=$3

INTERVAL=1
START_TIME=$(date +%s)

# Header only if file is new
if [ ! -f "$OUTFILE" ]; then
    echo "sample,timestamp,mem_MB" > "$OUTFILE"
fi

while kill -0 "$PID" 2>/dev/null; do
    PIDS=$(pstree -p "$PID" | grep -oP '\(\d+\)' | grep -oP '\d+')
    PIDS="$PID $PIDS"

    MEM_KB=$(ps -o rss= -p $PIDS 2>/dev/null | awk '{sum+=$1} END{print sum}')
    MEM_MB=$(echo "scale=2; $MEM_KB / 1024" | bc)

    NOW=$(date +%s)
    ELAPSED=$((NOW - START_TIME))

    echo "$SAMPLE,$ELAPSED,$MEM_MB" >> "$OUTFILE"

    sleep $INTERVAL
done
