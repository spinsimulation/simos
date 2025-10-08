#!/usr/bin/env bash
set -euo pipefail

PORT=${PORT:-8080}
BUILD_DIR="_build"
URL="http://localhost:${PORT}"

# Build docs
sphinx-build -v . "$BUILD_DIR"

# Open browser (portable) after a short delay, in the background
( sleep 1; python -m webbrowser "$URL" >/dev/null 2>&1 || true ) &

# Serve the site in the foreground (so Ctrl-C stops it)
python -m http.server "$PORT" -d "$BUILD_DIR"