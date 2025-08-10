.PHONY: fmt lint check test smoke bench doc clean

fmt:
\t@./scripts/rs_preflight.sh fmt

lint:
\t@./scripts/rs_preflight.sh lint

# check:
# \t@./scripts/rs_preflight.sh check

smoke:
\t@./scripts/rs_preflight.sh test --smoke

# test:
# \t@./scripts/rs_preflight.sh test --full

bench:
\t@cargo bench

doc:
\t@RUSTDOCFLAGS="-D warnings" cargo doc --no-deps

clean:
\t@cargo clean
