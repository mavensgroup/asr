#!/usr/bin/env bash
# ==============================================================================
# build_docs.sh — Build and serve asr documentation (mdBook + FORD)
#
# Usage:
#   ./tools/build_docs.sh           Build both mdBook manual and FORD Fortran API
#   ./tools/build_docs.sh --book    Build only the mdBook manual
#   ./tools/build_docs.sh --ford    Build only the FORD Fortran API
#   ./tools/build_docs.sh --serve   Serve mdBook locally with live-reload
#   ./tools/build_docs.sh --watch   Watch mdBook and rebuild on changes
# ==============================================================================

set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

MODE="${1:-all}"

build_book() {
    echo "==> [mdBook] Building manual from doc/manual/..."
    [ -e doc/manual/api ] || ln -sfn ../ford_api doc/manual/api
    mdbook build doc
    echo "==> [mdBook] Manual built at doc/build/book/index.html"
}

build_ford() {
    echo "==> [FORD] Parsing Fortran code in src/ and app/..."
    ford ford.md
    echo "==> [FORD] API reference built at doc/ford_api/index.html"
}

case "$MODE" in
    --book)
        build_book
        ;;
    --ford)
        build_ford
        ;;
    --serve)
        [ -d doc/ford_api ] || build_ford
        [ -e doc/manual/api ] || ln -sfn ../ford_api doc/manual/api
        echo "==> Serving documentation with live reload (http://localhost:3000)..."
        echo "==> Press Ctrl+C to stop."
        mdbook serve doc --open
        ;;
    --watch)
        [ -e doc/manual/api ] || ln -sfn ../ford_api doc/manual/api
        echo "==> Watching doc/manual/ for changes..."
        mdbook watch doc
        ;;
    all|"")
        build_ford
        build_book
        echo ""
        echo "================================================================="
        echo "Documentation build complete!"
        echo "  - User Manual:    file://$REPO_ROOT/doc/build/book/index.html"
        echo "  - Fortran API:    file://$REPO_ROOT/doc/build/book/api/index.html"
        echo "To browse locally with live-reload while writing:"
        echo "  ./tools/build_docs.sh --serve"
        echo "================================================================="
        ;;
    *)
        echo "Unknown option: $MODE"
        echo "Usage: $0 [--book | --ford | --serve | --watch]"
        exit 1
        ;;
esac
