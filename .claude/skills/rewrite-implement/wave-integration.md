# Cross-Compiling for Linux x86_64

The rewritten binary must ultimately run on Linux x86_64 (standard cloud executor target).
If developing on macOS or another platform, cross-compile using one of these options.

For full pipeline integration details (Nextflow process definition, Wave config, module
structure, testing), see the `/rewrite-integrate` skill.

## Option 1: Using `cross` (recommended)
```bash
cargo install cross
cross build --release --target x86_64-unknown-linux-gnu
```

## Option 2: Using Docker directly
```bash
docker run --rm -v "$PWD":/app -w /app rust:latest \
  cargo build --release --target x86_64-unknown-linux-gnu
```

## Option 3: Static linking with musl
```bash
rustup target add x86_64-unknown-linux-musl
cargo build --release --target x86_64-unknown-linux-musl
```
Note: `rust-htslib` requires additional setup for musl (static htslib build).

## Option 4: Zig as linker
```bash
cargo install cargo-zigbuild
cargo zigbuild --release --target x86_64-unknown-linux-gnu
```

## Verify the binary

Check the binary is statically linked or only depends on common system libraries:
```bash
# Inside a Linux environment
ldd ./target/x86_64-unknown-linux-gnu/release/<binary>
# Ideally shows only libc, libm, libpthread, libdl
```
