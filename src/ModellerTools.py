#!/usr/bin/python3

import sys
import os
import shutil

from modeller import *
from modeller.scripts import complete_pdb

def add_atoms(pdbfile):
    code    = pdbfile.split('.')[0]
    env                     = environ()
    env.edat.dynamic_sphere = True
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    mdl         = complete_pdb(env,pdbfile)
    patched_pdb = code+'_add.pdb'
    mdl.write(file=patched_pdb)

#############MAIN##############
if __name__ == "__main__":
    pdbfile = sys.argv[1]
    add_atoms(pdbfile)
