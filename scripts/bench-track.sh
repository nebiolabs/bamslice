#!/bin/bash
# Track benchmark results over time
set -e

HISTORY_DIR="benches/history"
CSV_FILE="$HISTORY_DIR/results.csv"
READABLE_FILE="$HISTORY_DIR/results.txt"
TEMP_FILE=$(mktemp)

# Get git info
COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
BRANCH=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")
TIMESTAMP=$(date "+%Y-%m-%d %H:%M:%S")
UNIX_TIME=$(date +%s)

# Run benchmarks and capture output
echo "Running benchmarks..."
cargo bench --bench bamslice_bench -- --output-format bencher 2>&1 | tee "$TEMP_FILE"

# Create CSV header if file doesn't exist
mkdir -p "$HISTORY_DIR"
if [ ! -f "$CSV_FILE" ]; then
    echo "timestamp,unix_time,commit,branch,benchmark,mean_ns,stddev_ns,throughput_mibs" > "$CSV_FILE"
fi

# Parse bencher format and extract results
# Format: "test process_blocks/small_500kb ... bench:     9225357 ns/iter (+/- 616480)"
grep "^test process_blocks" "$TEMP_FILE" | while read -r line; do
    # Extract benchmark name
    BENCH_NAME=$(echo "$line" | awk '{print $2}' | sed 's|process_blocks/||')

    # Extract mean time (ns/iter)
    MEAN_NS=$(echo "$line" | grep -o '[0-9]\+ ns/iter' | grep -o '[0-9]\+')

    # Extract stddev
    STDDEV_NS=$(echo "$line" | grep -o '+/- [0-9]\+' | grep -o '[0-9]\+')

    # Calculate throughput based on benchmark name
    case "$BENCH_NAME" in
        small_500kb)
            BYTES=500000
            ;;
        medium_1mb)
            BYTES=1000000
            ;;
        large_5mb)
            BYTES=5000000
            ;;
        *)
            BYTES=0
            ;;
    esac

    # Calculate MiB/s: (bytes / (1024*1024)) / (ns / 1e9) = bytes * 1e9 / (ns * 1048576)
    if [ "$MEAN_NS" -gt 0 ] && [ "$BYTES" -gt 0 ]; then
        THROUGHPUT=$(awk "BEGIN {printf \"%.2f\", ($BYTES * 1000000000.0) / ($MEAN_NS * 1048576.0)}")
    else
        THROUGHPUT="0.00"
    fi

    # Append to CSV
    echo "$TIMESTAMP,$UNIX_TIME,$COMMIT,$BRANCH,$BENCH_NAME,$MEAN_NS,$STDDEV_NS,$THROUGHPUT" >> "$CSV_FILE"
done

# Also maintain readable format
echo "" >> "$READABLE_FILE"
echo "=== $TIMESTAMP | $COMMIT | $BRANCH ===" >> "$READABLE_FILE"
grep "^test" "$TEMP_FILE" >> "$READABLE_FILE" || true
echo "----------------------------------------" >> "$READABLE_FILE"

# Clean up
rm -f "$TEMP_FILE"

echo ""
echo "Results appended to:"
echo "  CSV: $CSV_FILE"
echo "  Human-readable: $READABLE_FILE"
echo ""
echo "Recent results (last 5 runs):"
tail -16 "$CSV_FILE" | column -t -s,
