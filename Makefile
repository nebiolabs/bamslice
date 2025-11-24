.PHONY: clippy lint test bench

clippy:
	cargo clippy --all-targets --all-features -- -D warnings

lint: clippy

test:
	cargo test --all-features

bench:
	./scripts/bench-track.sh

flamegraph:
	./scripts/flamegraph.sh

coverage:
	./scripts/coverage.sh
