# -*- coding: utf-8 -*-
# conformational sampling analysis

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns

import tqdm

import pathlib
import os.path
import itertools

import alchemlyb.parsing.gmx




def read_xvg(filepath):
    """Parse XVG file `filepath` into a DataFrame."""
    df = alchemlyb.parsing.gmx._extract_dataframe(filepath)
    df.reset_index(inplace=True)
    df.rename(columns={df.columns[0]: "time", df.columns[1]: "angle"}, inplace=True)
    return df


def get_metadata(filename):
    """Extract metadata information from standardized file name.

    Parameters
    ----------
    `filename` : str

       must look like ::

            {forcefield}-{molid}-{solvent}-{interaction}-{lambda>}-{dih}-ts.xvg.bz2

    Returns
    -------
    tuple
           molid, forcefield, solvent, interaction, lambda, dihnumber


    Example
    -------
    Use to build a dataframe of all files::

       datafiles = pd.DataFrame([get_metadata(fn) for fn in
                                 itertools.chain(*(d.glob("*-ts.xvg*") for d in DATADIRS))],
                                 columns=["molid", "forcefield", "solvent", "interaction",
                                          "lambda", "dihedral", "filename"])
       datafiles = datafiles.sort_values(by=["molid", "forcefield", "dihedral", "solvent",
                                             "interaction", "lambda"]).reset_index(drop=True)

    """
    # specific for our simulations
    forcefields = ("cgenff", "oplsaa", "ligpargen", "gaff")
    solvents = ("water", "octanol")
    interactions = ("Coulomb", "VDW")


    # molid, forcefield, solvent, interaction, lambda, dihnumber
    filename = pathlib.Path(filename)
    forcefield, molid, solvent, interaction, lmbda, dih, _ = filename.name.split("-")
    for value, allowed in zip([forcefield, solvent, interaction], [forcefields, solvents, interactions]):
        if value not in allowed:
            raise ValueError(f"{filename}: {value} not in {allowed})")
    return molid, forcefield, solvent, interaction, float(lmbda)/1000, dih, filename



def periodic_angle(df, padding=45):
    """pad with ±padding beyond the edges and return angle series

    Parameters
    ----------
    df : DataFrame
       Must contain at least ``angle`` column (in degrees), in the
       range -180 .. 180 degrees.

    padding : float
       Extend the time series by +/- `padding` degrees by repeating
       data.

    Returns
    -------
    Series
       Sorted series of the angle with padding.

    """

    df = df.angle.sort_values()
    return pd.concat([df[df > 180 - padding] - 360, df, df[df < -180 + padding] + 360]).reset_index(drop=True)

def extract(datafiles, start=None, stop=None, step=None, padding=45):
    """extract angle data from GROMACS xvg file collection with periodic padding in long form

    Parameters
    ----------
    datafiles : DataFrame
       DataFrame with rows built from meta data with
       ``columns=["molid", "forcefield", "solvent", "interaction", "lambda", "dihedral", "filename"]``
    start : int or None
       start frame number (0-based)
    stop : int or None
       stop frame number (exclusive)
    step : int or None
       include every `step` frame
    padding : float
       periodically extend data by +/-`padding` degrees (see :func:`periodic`)

    Returns
    -------
    DataFrame
       Long-form DataFrame where each row corresponds to a single dihedral measurement
       in column *angle*;
       ``columns=["solvent", "interaction", "lambda", "dihedral", "angle"]``


    """
    dataframes = []
    # get each row as a df
    for i in tqdm.tqdm(range(len(datafiles))):
        simulation = datafiles[i:i+1]
        a = periodic_angle(read_xvg(simulation.iloc[0].filename).iloc[start:stop:step], padding=padding)
        identifiers = simulation.drop(columns=["molid", "forcefield", "filename"])
        df = pd.concat(
               [pd.concat(len(a)*[identifiers], axis="index")
                .reset_index(drop=True), a], axis="columns")
        dataframes.append(df)
    return pd.concat(dataframes)



def dihedral_violins(df, width=0.9):
    """Plot distributions of all dihedrals as violin plots.

    Parameters
    ----------
    df : DataFrame
         long form dataframe with lambda, solvent, and interaction
    width : float
          width of the violin element (>1 overlaps)

    Returns
    -------
    seaborn.FacetGrid
    """
    # number of Coul windows + 1 / number of VDW windows
    # (+1 for additional space with axes)
    width_ratios = [len(df[df['interaction'] == "Coulomb"]["lambda"].unique()) + 1,
                    len(df[df['interaction'] == "VDW"]["lambda"].unique())]

    g = sns.catplot(data=df, x="lambda", y="angle", hue="solvent", col="interaction",
                    kind="violin", split=True, width=width, inner=None, cut=0,
                    linewidth=0.5,
                    hue_order=["water", "octanol"], col_order=["Coulomb", "VDW"],
                    sharex=False, sharey=True,
                    height=2, aspect=2,
                    facet_kws={'ylim': (-180, 180),
                               'gridspec_kws': {'width_ratios': width_ratios,
                                                # 'wspace': 0.03
                                                }})
    g.set_xlabels(r"$\lambda$")
    g.set_ylabels(r"dihedral angle $\phi$")
    g.despine(offset=5)

    axC = g.axes_dict['Coulomb']
    axC.yaxis.set_major_locator(plt.matplotlib.ticker.MultipleLocator(60))
    axC.yaxis.set_minor_locator(plt.matplotlib.ticker.MultipleLocator(30))
    axC.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter(r"$%g^\circ$"))

    axV = g.axes_dict['VDW']
    axV.yaxis.set_visible(False)
    axV.spines["left"].set_visible(False)

    return g

def plot_violins(datafiles, molid, forcefield, dihedral,
                 start=None, stop=None, step=100,
                 figdir=None,
                 **kwargs):
    """Plot dihedral violins from `datafiles`.

    This function extracts the data and plots them with
    :func:`dihedral_violins`. It only needs as input the metadata
    DataFrame and a selection of `molid`, `forcefield` and the
    `dihedral` to plot.

    Parameters
    ----------
    datafiles : DataFrame
       List of input files with associated metadata.  `molid`,
       `forcefield`, and `dihedral` must be values in columns of the
       same name in the DataFrame

    molid : str
       Molecular ID to analyze, e.g. "SM46"

    forcefield : str
       Force field to analyze, e.g. "cgenff" or "oplsaa"

    dihedral : str
       Dihedral to analze, e.g. "dih1" or "dih2"

    start : int or None
       start frame number (0-based)

    stop : int or None
       stop frame number (exclusive)

    step : int or None
       include every `step` frame

    figdir : str or pathlike
       If not ``None``, write PDF plot to
       ``{figdir}/{molid}_{forcefield}_{dihedral}_violins.pdf``.

    kwargs : kwargs
        additional kwargs for :func:`dihedral_violins`

    Returns
    -------
    seaborn.FacetGrid

    """
    g = datafiles.groupby(by=["molid", "forcefield", "dihedral"])

    files = g.get_group((molid, forcefield, dihedral))

    df = extract(files, start=start, stop=stop, step=step)

    g = dihedral_violins(df, **kwargs)

    if figdir is not None:
        figfile = pathlib.Path(figdir) / f"{molid}_{forcefield}_{dihedral}_violins.pdf"
        g.savefig(figfile)
        #print(f"Violin plot saved to {figfile}.")

    return g
