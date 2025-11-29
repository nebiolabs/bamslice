.PHONY: clippy lint test bench perf coverage

clippy:
	cargo clippy --all-targets --all-features -- -D warnings

lint: clippy

test:
	cargo test --all-features

bench:
	./scripts/bench-track.sh

coverage:
	./scripts/coverage.sh

perf:
	cargo build --release && samply record ./target/release/bamslice -i tests/fixtures/dnbseq-test1.bam -s 0 -e 100000000 --log-level ERROR > /dev/null 2>/dev/null
