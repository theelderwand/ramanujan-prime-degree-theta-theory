#!/bin/bash
# Build the paper PDF. Requires pdflatex (texlive-latex-base + recommended).
# Run twice so cross-references resolve.
set -e
cd "$(dirname "$0")"
pdflatex -interaction=nonstopmode paper.tex
pdflatex -interaction=nonstopmode paper.tex
echo "Built: $(ls -la paper.pdf)"
echo "Overfull boxes: $(grep -c Overfull paper.log || true)"
echo "Undefined refs/cites: $(grep -c 'undefined' paper.log || true)"
