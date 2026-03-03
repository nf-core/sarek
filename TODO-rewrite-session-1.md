# Rewrite Session 1 (2026-03-03)

## Scope

- Phase 1 implementation: CLI scaffold + deterministic placeholder output plumbing for `ironqc`.

## Completed

- Created new Rust project structure under `ironqc/`.
- Added `Cargo.toml` with requested dependencies:
  - `clap` (derive)
  - `anyhow`
  - `rust-htslib` with `static`, `bzip2`, `lzma`, `curl`
  - `indexmap`
  - `rayon`
  - `plotters`
- Added release profile tuning:
  - `lto = true`
  - `codegen-units = 1`
  - `strip = true`
  - `opt-level = 3`
- Implemented CLI scaffold with subcommands:
  - `stats`
  - `mosdepth`
  - `indexcov`
  - `bundle`
- Implemented placeholder output file creation modules:
  - `src/stats/mod.rs`
  - `src/mosdepth/mod.rs`
  - `src/indexcov/mod.rs`
  - `src/bundle/mod.rs`
  - shared filesystem helpers in `src/io.rs`
- Added integration test skeleton in `ironqc/tests/integration_cli.rs` to validate expected files are created.

## Pipeline Invocation Verification Performed

- Checked `modules/nf-core/samtools/stats/main.nf` for stats output naming.
- Checked `modules/nf-core/mosdepth/main.nf` for argument semantics (`--by`, `-n`, `--fast-mode`).
- Checked `modules/nf-core/goleft/indexcov/main.nf` and `conf/modules/indexcov.config` for indexcov output placement and naming.

## Blockers

- `cargo` is not installed in this environment (`command not found`), so these required checks could not be executed:
  - `cargo fmt`
  - `cargo clippy -- -D warnings`
  - `cargo test`
- `rust-analyzer` is not installed, so `lsp_diagnostics` could not run for changed Rust files.

## Next Steps

1. Install Rust toolchain (`cargo`, `rustc`) and `rust-analyzer`.
2. Run `cargo fmt && cargo clippy -- -D warnings && cargo test` in `ironqc/`.
3. Fix any resulting issues and re-run until clean.
4. Continue Phase 1 acceptance validation against invocation expectations.

## Remaining Estimate

- Phase 1 remaining effort: ~1 session (toolchain setup + verification pass + final cleanup).
