"""
simulate.py

Performs 

Updated to python3:
-
"""

import os
import sys
import re

import OsTools as ot
import MFFTools as mfft
import PullSim as ps

def do_simulations(pdbfile,	# input PDB file
                   pair,	# set of chains
                   distances,	# list of distances
                   boRun=True,	# start production run after equilibration
                   pull='cons', # type of pull: _cons_traints or _umbrella_
                   nsteps=None  # nsteps for production run (None=default)
                   ):

    # make directory for output:
    pdbase  = pdbfile.split('.pdb')[0]
    dir     = pdbase
    if (not (os.path.exists(dir) and os.path.isdir(dir))):
        os.mkdir(dir)
    os.chdir(dir)
    # create symlink to input pdbfile
    ot.symlink('../'+pdbfile,pdbfile)

    # create CG (coarse-grained MARTINI) representation of protein:
    cggro_merged, cgndx_merged, cg_stable_itppair, \
        gro_cgpair, nat_cgpair, fgpair = \
            mfft.do_fg2cg(pdbfile, pair)

    print("After CG2FG got output files:", cggro_merged, cgndx_merged, 
          cg_stable_itppair, gro_cgpair, nat_cgpair, fgpair)

    groups = []
    for i in gro_cgpair:
        groups.append(re.sub('_cg.gro','',i))

    root_dir = os.getcwd()
    for distance in sorted(distances):
        result = ps.do_pullsim(root_dir,distance,cggro_merged,cgndx_merged,
                               groups,cg_stable_itppair,
                               nat_cgpair,pdbase,fgpair,gro_cgpair,
                               boRun,pull,nsteps)
        if result:
            print("pullsim failed in", os.getcwd())
            # make empty file to flag this error
            open("FAILED", 'w').close()

    os.chdir(root_dir)

#last line
