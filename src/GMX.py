"""
GMX.py

Contains functions that call specific gromacs commands with options/flags.

Modified for python 3
"""
import os
import re

import OsTools as ot
import Env as env

def g_tpbconv(*options, **keywords):
    return ot.system(env.g_tpbconv, *options, **keywords)

def g_gmxdump(*options, **keywords):
    return ot.system(env.g_gmxdump, *options, **keywords)

def g_energy(*options, **keywords):
    return ot.system(env.g_energy, *options, **keywords)

def g_pdb2gmx(*options, **keywords):
    return ot.system(env.pdb2gmx, *options, **keywords)

def g_editconf(*options, **keywords):
    return ot.system(env.editconf, *options, **keywords)

def g_make_ndx(*options, **keywords):
    return ot.system(env.make_ndx, *options, **keywords)

def gen_ndxUIlist(offset,gro_cgpair,nat_cgpair):
    stdinlist = []
    count = 0
    for i in range(len(gro_cgpair)):
        start = str(count+1)
        count = nat_cgpair[i]+int(count)
        end   = str(count)
        stdinlist.append('\"Protein\" & a '+start+'-'+end+'')
        group = re.sub('_cg.gro','',gro_cgpair[i])
        stdinlist.append('name '+str(offset)+' '+group+'')
        offset += 1
        stdinlist.append('\"'+group+'\" & a B*'+'')
        group += '_backbone'
        stdinlist.append('name '+str(offset)+' '+group+'')
        offset += 1
    stdinlist.append('q')
    return stdinlist

def g_g_mindist(*options, **keywords):
    return ot.system(env.g_mindist, *options, **keywords)

def g_genrestr(*options, **keywords):
    return ot.system(env.genrestr, *options, **keywords)

def g_genbox(*options, **keywords):
    return ot.system(env.genbox, *options, **keywords)

def g_grompp(*options, **keywords):
    return ot.system(env.grompp, *options, **keywords)

def g_mdrun(*options, **keywords):
    return ot.system(env.mdrun, *options, **keywords)

def g_analyze(*options, **keywords):
    return ot.system(env.g_analyze, *options, **keywords)
