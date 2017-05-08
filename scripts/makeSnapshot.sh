#!/bin/bash
set -euo pipefail
NAME=CAW-$(git describe --tags)
git archive HEAD --prefix=$NAME/ | gzip > $NAME.tar.gz
echo "Wrote $NAME.tar.gz"
