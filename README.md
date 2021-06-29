# Conformational Analysis of SAMPL7 molecules

Start with SM46 (see submitted paper).

Bogdan created timeseries for the two main dihedrals (see tar file on spudda,
too big for repo).

To set up:

    bash ./bin/RSYNC
    bash ./bin/SETUP

You can then use the notebooks in analysis.

You need a conda environment with

- alchemlyb
- pandas
- matplotlib
- numpy
- tqdm
- seaborn



## Raw data

Generated with `gmx angle` (see `GROMACS_ANALYZE_dihedral_analysis`):
```bash
  gmx angle -type dihedral -f $i/md.xtc -n $ff-dih1.ndx -od $j-dih1-dist.xvg -ov $j-dih1-ts.xvg
  gmx angle -type dihedral -f $i/md.xtc -n $ff-dih2.ndx -od $j-dih2-dist.xvg -ov $j-dih2-ts.xvg
```

### CGENFF and OPLS-AA (mol2ff)

Tar file with data. 
* `spudda.beckstein.physics.asu.edu:/nfs/homes/bogdan/Projects/SAMPL7-PhysProp/_dih_analysis.tgz`

Use `bin/RSYNC.sh` to get the data.

### GAFF and LigParGen

Shujie generated the data

* `yamsolo.beckstein.physics.asu.edu:/nfs/homes3/sfan/Projects/SAMPL7/_dih_gaff`
* `yamsolo.beckstein.physics.asu.edu:/nfs/homes3/sfan/Projects/SAMPL7/_dih_ligpargen`
