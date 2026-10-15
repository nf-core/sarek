#!/usr/bin/env bash
set -euo pipefail

root=$(dirname "$(realpath "$0")")
revision=21b5ed2114207b93439e283d2bc2107e1420b4bf
checksum=e1f59765101586b08f2033d90c592553e8c3bfbb97b3ae0976b8921dd298db06
compiler=${CC:-cc}
command -v "$compiler" >/dev/null || { echo "Nextflow parser setup requires a C compiler (CC or cc)." >&2; exit 1; }

case "$(uname -s)-$(uname -m)" in
    Darwin-arm64|Darwin-x86_64) extension=dylib; link_flag=-dynamiclib ;;
    Linux-aarch64|Linux-x86_64) extension=so; link_flag=-shared ;;
    *) echo "Nextflow parser setup supports macOS/Linux on arm64/x86_64." >&2; exit 1 ;;
esac

temporary=$(mktemp -d)
trap 'rm -rf "$temporary"' EXIT
curl --fail --location --retry 3 \
    "https://codeload.github.com/nextflow-io/tree-sitter-nextflow/tar.gz/$revision" \
    -o "$temporary/parser.tar.gz"
if command -v sha256sum >/dev/null; then
    actual=$(sha256sum "$temporary/parser.tar.gz")
else
    actual=$(shasum -a 256 "$temporary/parser.tar.gz")
fi
if [[ ${actual%% *} != "$checksum" ]]; then
    echo "Nextflow parser source checksum mismatch." >&2
    exit 1
fi
tar -xzf "$temporary/parser.tar.gz" -C "$temporary"
source_dir="$temporary/tree-sitter-nextflow-$revision"
"$compiler" -O2 -std=c11 -fPIC "$link_flag" -I"$source_dir/src" \
    "$source_dir/src/parser.c" "$source_dir/src/scanner.c" \
    -o "$temporary/libnextflow.$extension"
mkdir -p "$root/lib"
cp "$source_dir/LICENSE" "$root/lib/LICENSE"
mv "$temporary/libnextflow.$extension" "$root/lib/libnextflow.$extension"
