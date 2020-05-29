#!/usr/bin/env python3
#PDB.py


"""
Author: Rene Pool (rpool@few.vu.nl)

Revision History:
File created on 24/11/2009

Parses a PDB file

Modified for python 3
"""

USAGE = """
Usage:
python PDB_ATOMS2FASTA.py ????.pdb [-splitpdb 0/1]
-----------------------------------------------------------------
Case -splitpdb 0: do not split pdb into separate chains (DEFAULT)
Case -splitpdb 1: split pdb into separate chains
-----------------------------------------------------------------
Output is in file ????.fasta or into separate fasta files for
split pdbs. If -splitpdb 1, we obtain the speparate PDB files
as well.
"""

import sys
import os
import string
from Bio.PDB import PDBParser,PPBuilder
import Dice_RPL # This own class overrides original residue numbering per chain
import PDBIO_RPL

seq3to1 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
           'CYX': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
           'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
           'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y',
           'VAL': 'V', 'ASX': 'B', 'GLX': 'Z'}

#def measure_peptide_bonds(infile):
#    parser = PDBParser()
#    struct = parser.get_structure('mystruct',infile)
#    ppb    = PPBuilder()
#    for pp in ppb.build_peptides(struct):
#        print(pp
#        print(pp.get_sequence()
#
##    for pp in ppb:
#
#
#    for model in struct:
#        for chain in model:
#            for index, pp in enumerate(ppb.build_peptides(chain)):
#                print(chain.id,index,pp.get_sequence(),pp
#
##            for resid in chain:
##                resid_it = iter(resid)
###                prev = resid_it.prev()
##                next = resid_it.next()
##                print(chain.id, next


#############RUN###############
def run(infile,splitpdb):
    parser = PDBParser()
    struct = parser.get_structure('mystruct',infile)
    ppb    = PPBuilder()

    basename = os.path.basename(infile)
    prefix = os.path.splitext(basename)[0]
    if splitpdb==0: # We do NOT split the PDB and fasta files!
        seqfile = open(prefix+'.fasta', 'w')
        pdbio=PDBIO_RPL.PDBIO()
        pdbio.set_structure(struct)
        cleanfile = prefix+'_clean.pdb'
        pdbio.save(cleanfile)
    ListChains = []
    for model in struct:
        for chain in model:
            ListChains.append(chain.id)
            ListPpdb = ppb.build_peptides(chain)
            if(len(ListPpdb)>0):
                for index, pp in enumerate(ListPpdb):
#                    print(chain.id,index,pp.get_sequence(),pp
                    if splitpdb==1: # We split the PDB and fasta files!
                        seqfile = open(prefix+'_'+chain.id+'.'+str(index)+'.fasta', 'w')
                    seq = pp.get_sequence()
                    seqfile.write('>%s %s\n' % (prefix+'_chain_'+chain.id+'_'+str(index), len(seq)))
                    seqfile.write('%s' % seq)
                    seqfile.write('\n')
                    if splitpdb==1: # We split the PDB and fasta files!
                        seqfile.close()
                        startres = pp[0].id[1]
                        endres   = pp[-1].id[1]
                        ofile = prefix+'_'+chain.id+'.'+str(index)+'.pdb'
                        Dice_RPL.extract(struct,chain.id,startres,endres,ofile)
            else:
#               Also split chains that do not consist of amino acids!
                ChainList = chain.get_list()
                startres  = ChainList[0].id[1]
                endres    = ChainList[0].id[-1]
                ofile = prefix+'_'+chain.id+'.'+str(index)+'.pdb'
                Dice_RPL.extract(struct,chain.id,startres,endres,ofile)
    if splitpdb==0: # We do NOT split the PDB and fasta files!
        seqfile.close()

    return ListChains

#    seq = []
#    for model in struct:
#        nchn = 0
#        for chain in model:
#            seq.append([])
#            for residue in chain:
##                skip HETATMS!
#                if residue.id[0]==' ':
#                    aa = residue.get_resname()
#                    seq[nchn].append(seq3to1.get(aa, aa))
#            nchn += 1

#    basename = os.path.basename(infile)
#    prefix = os.path.splitext(basename)[0]
#    if splitpdb==0: # We do NOT split the PDB and fasta files!
#        seqfile = open(prefix+'.fasta', 'w')
#    for model in struct:
#        nchn = 0
#        for chain in model:
#            if splitpdb==1: # We split the PDB and fasta files!
#                seqfile = open(prefix+'_'+chain.id+'.fasta', 'w')
#            seqfile.write('>%s %s\n' % (prefix+'_chain_'+chain.id, len(seq[nchn])))
#            seq[nchn] = string.join(seq[nchn], '')
#            seq[nchn] = string.strip(seq[nchn])
#            seqfile.write('%s' % seq[nchn])
#            seqfile.write('\n')
#            nchn += 1
#            if splitpdb==1: # We split the PDB and fasta files!
#                seqfile.close()
#                l = []
#                for residue in chain:
#                    if residue.id[0]==' ':
#                        l.append(residue.id[1])
#                if len(l)>1:
#                    start = min(l)
#                    end   = max(l)
#                    ofile = prefix+'_'+chain.id+'.pdb'
#                    extract(struct,chain.id,start,end,ofile)
#
#    if splitpdb==0: # We do NOT split the PDB and fasta files!
#        seqfile.close()

#############MAIN##############
if __name__ == "__main__":

    splitpdb = 0 # TAG for optional splitting od PDB file and FASTA files
    narg = len(sys.argv)
    if narg<2:
        print(USAGE)
        sys.exit()
    elif narg<3:
        splitpdb = 0
    elif narg<4:
        print("Wrong number of input parameters!")
        print(USAGE)
        sys.exit()
    elif narg==4:
        if sys.argv[2]=='-splitpdb':
            tmparg = int(sys.argv[3])
            if (tmparg==0 or tmparg==1):
                splitpdb = tmparg
            else:
                print("Wrong input for -splitpdb")
                print(USAGE)
                sys.exit()
        else:
            print("Wrong syntax of input parameters!")
            print(USAGE)
            sys.exit()

    infile = sys.argv[1]

#    measure_peptide_bonds(infile)

    run(infile,splitpdb)



