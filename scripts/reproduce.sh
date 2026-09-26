#!/usr/bin/env sh
# Reproduce every product of the paper companion (tests -> figures --no-cache -> paper numbers
# -> provenance check -> acceptance [-> tectonic if available]).  Thin wrapper around
# scripts/reproduce.py, which works the same on Windows, Linux and macOS.
#
# Usage:  sh scripts/reproduce.sh [--quick] [--with-sweep] [--skip-tests] [--no-latex] [--tectonic PATH]
set -e
cd "$(dirname "$0")/.."
# $PYTHON wins; otherwise python3, unless it is missing or cannot run (on Windows `python3` can be
# the Microsoft Store stub, which exists on PATH but does not start Python), then python.
if [ -n "${PYTHON:-}" ]; then
    PY="$PYTHON"
elif command -v python3 >/dev/null 2>&1 && python3 -c "import sys" >/dev/null 2>&1; then
    PY=python3
else
    PY=python
fi
exec "$PY" scripts/reproduce.py "$@"
