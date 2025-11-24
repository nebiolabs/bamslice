sudo CARGO_PROFILE_RELEASE_DEBUG=true RUST_LOG=info cargo flamegraph -- -i tests/fixtures/dnbseq-test1.bam -s 0 -e 100000000 --log-level ERROR > /dev/null 2>/dev/null
