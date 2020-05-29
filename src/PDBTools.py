"""
PDBTools.py

Contains different tools for working with pdb formatted files.


Modified for python3
"""
import os
import re

import Env as env
import OsTools as ot

def merge_pdbs(base,pair):
    print("Starting merge:", base, pair)
    ls  = os.listdir(os.getcwd())
    pat = re.compile(base+'_.\.[0-9].*\.pdb') # e.g. '1VET_A.0.pdb'
    pdbls = []
    for dirfile in ls:
        if pat.match(dirfile):
            pdbls.append(dirfile)
    print("Merging", sorted(pdbls))
    # we're supposed to get:
    # fg1.pdb fg1_A.0.pdb fg1_A.0_cg.pdb fg1_B.0.pdb fg1_B.0_cg.pdb
    # but second time in same dir, we run into our own output files:
    # fg1_A.pdb fg1_A_cg.pdb fg1_B.pdb fg1_B_cg.pdb
    pairfiles = []
    for chain in pair:
        chainsplit = chain.split(',')
        fname  = base+'_'+''.join(chainsplit)+'.pdb'
        f      = open(fname,'w')
        pairfiles.append(fname)
        for chain in chainsplit:
            p = re.compile('_'+chain+'\.')
            for file in sorted(pdbls):
                if (p.search(file)):
                    ff = open(file,'r')
                    f.write(ff.read())
                    ff.close()
        f.close()
    return pairfiles

def merge_fastas(base,pair):
    ls      = os.listdir(os.getcwd())
    ext     = re.compile('fasta')
    fastals = []
    for i in ls:
        if (ext.match(i.split('.')[-1])):
            fastals.append(i)
    pairfiles = []
    for i in pair:
        isplit = i.split(',')
        fname  = base+'_'+''.join(isplit)+'.fasta'
        f      = open(fname,'w')
        pairfiles.append(fname)
        for j in isplit:
            p = re.compile('_'+j+'\.')
            for file in sorted(fastals):
                if (p.search(file)):
                    ff = open(file,'r')
                    f.write(ff.read())
                    ff.close()
        f.close()
    return pairfiles

def pdb2dssp(fgpair):
    dssppair = []
    for fgfile in fgpair:
        dsspname = re.sub('.pdb','.dssp',fgfile)
        if ot.system(env.dssp, options=[fgfile, dsspname], log='dsspcmbi.err'):
            print("dssp failed for", fgfile)
            exit(-1)
        dssppair.append(dsspname)
    return dssppair

def check_missing_atoms(pdbfile):
    boMissingAtoms = False
    fr = open(pdbfile,'r')
    pdb_mem = fr.readlines()
    fr.close()
    for line in pdb_mem:
        if(re.search('MISSING',line) or re.search('missing',line)):
            if(re.search('ATOM',line) or re.search('atom',line)):
                boMissingAtoms = True
        if(re.search('REMARK 470',line)):
            boMissingAtoms = True
    del pdb_mem
    return boMissingAtoms
