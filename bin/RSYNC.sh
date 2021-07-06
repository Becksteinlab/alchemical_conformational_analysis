rsync -avP spudda.beckstein.physics.asu.edu:/nfs/homes/bogdan/Projects/SAMPL7-PhysProp/_dih_analysis.tgz .

rsync -avP --exclude="*-dist.xvg.bz2" yamsolo.beckstein.physics.asu.edu:/nfs/homes3/sfan/Projects/SAMPL7/_dih_gaff .

rsync -avP --exclude="*-dist.xvg.bz2" yamsolo.beckstein.physics.asu.edu:/nfs/homes3/sfan/Projects/SAMPL7/_dih_ligpargen .

# rename these files in bin/SETUP.py
rsync -avP spudda.beckstein.physics.asu.edu:/nfs/homes/bogdan/Projects/SAMPL7-PhysProp/dih*ts_SM46_water_VDW_1000_extended_1microsec.xvg.gz .
rsync -avP spudda.beckstein.physics.asu.edu:/nfs/homes/bogdan/Projects/SAMPL7-PhysProp/dih[12]ts_SM46_octanol_VDW_1000_extended_466ns.xvg.gz .
