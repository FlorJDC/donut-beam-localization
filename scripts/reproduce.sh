#!/usr/bin/env sh
# Reproduce every product of the paper companion (tests -> figures --no-cache -> paper numbers
# -> provenance check -> acceptance [-> tectonic if available]).  Thin wrapper around
# scripts/reproduce.py, which works the same on Windows, Linux and macOS.
#
# Usage:  sh scripts/reproduce.sh [--quick] [--with-sweep] [--skip-tests] [--no-latex] [--tectonic PATH]
set -e
cd "$(dirname "$0")/.."
PY="${PYTHON:-python3}"
command -v "$PY" >/dev/null 2>&1 || PY=python
exec "$PY" scripts/reproduce.py "$@"
