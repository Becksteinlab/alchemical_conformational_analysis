#!/bin/bash

TARFILE=_dih_analysis.tgz
PARALLEL=$(dirname $0)/parallel.py

tar xvf $TARFILE || { echo "Missing tar file $TARFILE"; exit 1; }

echo "removing histograms (not needed)"
find . -name '*-dist.xvg' | xargs rm

echo "compressing xvgs to save space (using 4 cores)"
find . -name '*.xvg' | xargs ${PARALLEL} 4 bzip2 -v ---
