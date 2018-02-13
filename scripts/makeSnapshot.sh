#!/bin/bash
set -euo pipefail
NAME=Sarek-$(git describe --tags)
git archive HEAD --prefix=$NAME/ | gzip > $NAME.tar.gz
echo "Wrote $NAME.tar.gz"
