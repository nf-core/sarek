# Rust Project Scaffold Reference

## Directory Structure

```
<project-name>/
├── Cargo.toml
├── README.md
├── AGENTS.md                    # Instructions for AI agents working on this code
├── src/
│   ├── main.rs                  # Entry point, CLI dispatch
│   ├── cli.rs                   # clap derive argument parsing
│   ├── config.rs                # YAML/TOML config (if needed)
│   ├── io.rs                    # Shared I/O utilities
│   ├── <domain>/                # Domain module (e.g., rna/, variant/, align/)
│   │   ├── mod.rs               # Re-exports sub-modules
│   │   ├── <tool1>/             # Per-tool sub-module
│   │   │   ├── mod.rs
│   │   │   └── ...
│   │   └── <tool2>/
│   │       └── ...
│   └── <annotation>.rs          # Annotation parser (GTF, BED, GFF)
├── tests/
│   ├── integration_test.rs      # Integration tests vs reference outputs
│   ├── data/                    # Test input files (small)
│   └── expected/                # Reference outputs from upstream tools
├── benchmark/
│   ├── README.md                # Benchmark methodology and results
│   ├── input/                   # Benchmark input data (or symlinks)
│   └── run_benchmark.sh         # Reproducible benchmark script
├── scripts/
│   └── validate_*.py            # Output comparison scripts
└── .github/
    └── workflows/
        └── ci.yml               # CI: fmt, clippy, test
```

## Cargo.toml Template

```toml
[package]
name = "<binary-name>"
version = "0.1.0"
edition = "2021"
description = "<one-line description>"

[[bin]]
name = "<binary-name>"
path = "src/main.rs"

[dependencies]
# CLI
clap = { version = "4", features = ["derive"] }

# Error handling
anyhow = "1"

# Logging
log = "0.4"
env_logger = "0.11"

# Serialization (for config files)
serde = { version = "1", features = ["derive"] }
serde_yaml = "0.9"

# BAM/CRAM I/O (bioinformatics)
# rust-htslib = { version = "0.47", features = ["static"] }

# Parallelism
rayon = "1"

# Ordered maps (when insertion order matters)
# indexmap = "2"

# Interval trees (for genomic interval queries)
# coitrees = "0.4"

# Plotting
# plotters = { version = "0.3", default-features = false, features = ["svg_backend", "bitmap_backend", "line_series"] }

# Compression
# flate2 = "1"

# CSV output
# csv = "1"

[profile.release]
opt-level = 3
lto = true
codegen-units = 1
strip = true
```

## CI Template (.github/workflows/ci.yml)

```yaml
name: CI
on:
  push:
    branches: [main]
  pull_request:

jobs:
  test:
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest]
    runs-on: ${{ matrix.os }}
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - uses: Swatinem/rust-cache@v2
      - run: cargo test --release

  format:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - run: cargo fmt --check

  clippy:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - uses: Swatinem/rust-cache@v2
      - run: cargo clippy -- -D warnings
```

## Code Conventions

- **Formatting**: Default `rustfmt` (no config file). 4-space indent.
- **Errors**: `anyhow::Result<T>` everywhere. Propagate with `?`. Context with `.context()`.
- **No unwrap**: `unwrap()` / `expect()` only in tests or with safety comment.
- **Docs**: `//!` module docs on every file. `///` on all public items.
- **Derives**: `#[derive(Debug)]` on all structs. Add `Clone`, `Default` as needed.
- **Naming**: Types=CamelCase, functions=snake_case, constants=SCREAMING_SNAKE_CASE.
- **Tests**: Co-located `#[cfg(test)] mod tests { use super::*; ... }` in each file.
  Integration tests in `tests/` directory.

## AGENTS.md Template

Write an AGENTS.md that documents:
1. What the project does (2-3 sentences)
2. Build/test/lint commands
3. Project structure (directory tree with descriptions)
4. Code style conventions
5. Key dependencies and their purpose
6. CI pipeline
7. Notes for agents (gotchas, important invariants)
