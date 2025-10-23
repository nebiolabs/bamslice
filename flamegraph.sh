sudo CARGO_PROFILE_RELEASE_DEBUG=true RUST_LOG=info cargo flamegraph -- -i tests/fixtures/emseq-test1.bam -s 0 -e 10000000 --log-level ERROR 2>/dev/null
