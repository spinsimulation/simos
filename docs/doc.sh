#!/usr/bin/env bash
sphinx-build -v . _build
open "http://localhost:8080"
python -m http.server 8080 -d _build
