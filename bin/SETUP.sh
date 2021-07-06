#!/bin/bash

TARFILE=_dih_analysis.tgz
PARALLEL=$(dirname $0)/parallel.py

# oplsaa and cgenff
test -e $TARFILE && tar xvf $TARFILE || { echo "Missing tar file $TARFILE - SKIPPING"; }


echo "removing histograms (not needed)"
find . -name '*-dist.xvg' | xargs rm

echo "compressing xvgs to save space (using 4 cores)"
find . -name '*.xvg' | xargs ${PARALLEL} 4 bzip2 -v ---

# purge xvg
find . -name '*.xvg' | xargs rm

# extended sims
XDIR=_dih_oplsaa_extended
echo "extended windows --> $XDIR"

mkdir -p $XDIR
for dih in 1 2; do
    mv dih${dih}ts_SM46_water_VDW_1000_extended_1microsec.xvg.gz $XDIR/oplsaa-SM46-water-VDW-1000-dih${dih}-ts.xvg.gz
    mv dih${dih}ts_SM46_octanol_VDW_1000_extended_466ns.xvg.gz $XDIR/oplsaa-SM46-octanol-VDW-1000-dih${dih}-ts.xvg.gz
done
