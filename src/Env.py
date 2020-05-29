"""
Env.py

Defines a set of environment variables.

Modified for python 3
"""
import os
import os.path
import sys
import re

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

try:
    mffhome   = os.environ['MFFHOME']
except KeyError:
    dir  = os.path.abspath(sys.argv[0])
    file = 'MFF'
    while not os.path.isdir(dir+'/'+file) and len(dir)>1:
        # trim last dir:
        dir = os.path.split(dir)[0] # returns ['/path', 'file']
    mffhome = dir+'/'+file
    if not os.path.isdir(mffhome):
        print('Cannot find', file, 'up from', sys.argv[0])
        exit(-1)

mfflib    = mffhome+'/lib'
sys.path.append(mfflib)
mffFF     = mffhome+'/FF'
mffitp    = 'martini_v2.1.itp'
mffions   = 'martini_v2.0_ions.itp'
mffSO4    = 'martini_sulfate.itp'
mffFSC    = 'martini_fusicoccin.itp'
mffstruct = mffhome+'/structures'
mffmass   = 72.0 # amu
mffepsW   = 5.0 # kJ/mol
mffsigW   = 0.47 # nm
home      = os.environ['HOME']
pywork    = 'PyWork/trunk'
pysrc     = 'src'
#usrbin    = '/usr/bin'
usrbin    = os.environ['HOME']+'/bin'
#gmxbin    = home+'/progs/gromacs-4.0.5/bin/' # include trailing slash!
#gmxbin    = home+'/progs/GMX_4.0.5_mpi/bin/' # include trailing slash!
gmxbin    = os.path.dirname(which('mdrun')) + '/' # include trailing slash!

# GMX tools
pdb2gmx   = gmxbin+'pdb2gmx'
editconf  = gmxbin+'editconf'
make_ndx  = gmxbin+'make_ndx'
genrestr  = gmxbin+'genrestr'
genbox    = gmxbin+'genbox'
grompp    = gmxbin+'grompp'
mdrun     = gmxbin+'mdrun'
g_mindist = gmxbin+'g_mindist'
g_analyze = gmxbin+'g_analyze'
g_energy  = gmxbin+'g_energy'
g_gmxdump = gmxbin+'gmxdump'
g_tpbconv = gmxbin+'tpbconv'

# MFF tools
fg2cg     = mfflib+'/atom2cg.awk'
dssp      = mfflib+'/dsspcmbi'
seq2itp   = mfflib+'/seq2itp.pl'
sel_restr = mfflib+'/select_restr.pl'
import dssp2ssd_NO_CALC_CHAIN_BREAKS as to_ssd
import PDB_ATOMS2FASTA as a2f
