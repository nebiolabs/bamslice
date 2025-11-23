#!/bin/bash
# Code coverage script using grcov

set -e

# Clean previous coverage data
echo "Cleaning previous coverage data..."
find . -name "*.profraw" -delete
rm -rf target/coverage/
mkdir -p target/coverage

# Set environment variables for coverage
export CARGO_INCREMENTAL=0
export RUSTFLAGS="-Cinstrument-coverage"
export LLVM_PROFILE_FILE="target/coverage/annotate-reads-%p-%m.profraw"

echo "Building and running tests with coverage instrumentation..."
cargo build
cargo test

echo "Generating coverage report with grcov..."
grcov target/coverage \
    --binary-path target/debug/deps/ \
    --source-dir . \
    --output-types html,lcov \
    --branch \
    --ignore-not-existing \
    --ignore "/*" \
    --ignore "target/*" \
    --ignore "tests/*" \
    --excl-line "GRCOV_EXCL_LINE" \
    --excl-start "GRCOV_EXCL_START" \
    --excl-stop "GRCOV_EXCL_STOP" \
    --output-path target/coverage/

echo "Coverage report generated:"
echo "  HTML: target/coverage/html/index.html"
echo "  LCOV: target/coverage/lcov"

# Calculate coverage percentage
if command -v lcov &> /dev/null; then
    echo ""
    echo "Coverage Summary:"
    lcov --summary target/coverage/lcov 2>/dev/null | grep -E "(lines|functions|branches)" || echo "Install lcov for detailed summary"
fi