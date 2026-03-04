FROM rust:1.83-slim-bookworm AS builder
RUN apt-get update && apt-get install -y --no-install-recommends git ca-certificates && rm -rf /var/lib/apt/lists/*
WORKDIR /build
RUN git clone --depth 1 -b feature/rewrite-tools https://github.com/evanfloden/iron-sarek.git .
WORKDIR /build/ironqc
RUN cargo build --release

FROM debian:bookworm-slim
COPY --from=builder /build/ironqc/target/release/ironqc /usr/local/bin/ironqc
