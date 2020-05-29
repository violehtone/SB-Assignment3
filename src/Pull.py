#!/usr/bin/env python3
#Pull.py
"""
Author: Rene Pool (r.pool@few.vu.nl)

Revision History:
File created on 24/11/2009
Modified for python 3 on 22/01/2020

Displaces the atoms of two groups along the vector connecting their COMs
"""

USAGE = """
Usage:
python Pull.py -f [pdb_file] -c [chains] -d [<start>:<stop>:<increment>]
"""

import os
import re
import sys

import Simulate as sim
import PDBTools as pdbt
import ModellerTools as mt
import math
from optparse import OptionParser#, OptionGroup

######## COMMAND LINE / INPUT STUFF ##########

def parse_commandline():
    usage = "%prog [<pdbfile>] [<chains>] [<distances] [options]"
    description = \
        "%prog runs pre-equilibration and per-distance equilibration "\
        "for PMF calculation. It can optionally use fine-grained (atomistic) "\
        "or coarse-grained forcefield. It can optionally run production "\
        "simulations automatically. It can use either constraints or "\
        "umbrella sampling."
    parser = OptionParser(usage=usage, description=description)

    # pdbfile
    # pair
    # drangefile
    # -no_prod
    # -umbrella
    # -submit --> leave out
    # -nsteps
    # -fg
    parser.add_option("-f", "--pdb",  dest="pdbfile", metavar="<file>",
                      help="input PDB file")
    parser.set_defaults(pdbfile=None)
    parser.add_option("-c", "--chains",  dest="pair", metavar="<chIDs:chIDs>",
                      help="colon separated pair of set of chains, e.g. 'A:B' or 'AB:DEF' (%default)")
    parser.set_defaults(pair="A:B")
    
    parser.add_option("-d", "--distances",  dest="drangefile",
                      metavar="<file/range>",
                      help="distances to calculate, can be either the name of "\
                      "a file containing a list of explicit distances, or a "\
                      "colon separated range of (start:stop:step), e.g. "\
                      "'2.3:5.2:0.1' (%default)")
    parser.set_defaults(drangefile=None)
    
    parser.add_option("-p", "--pull", dest="pull", metavar="<type>",
                      help="Type of pulling to perform: cons for constraints, "\
                      "or umbr for umbrella sampling (%default)")
    parser.set_defaults(pull='umbr')
    
    parser.add_option("-r", "--production", dest="prod",
                      action="store_true",
                      help="Perform production runs automatically after "\
                      "equilibriation (%default)")
    parser.set_defaults(prod=False)
    
    parser.add_option("-i", "--nsteps", dest="nsteps", type="int", metavar="i",
                      help="Total number of steps for production simulations (%default)")
    parser.set_defaults(nsteps=None)
    
    parser.add_option("", "--fg", dest="fg",
                      action="store_true",
                      help="Use atomistic (fine-grained) forcefield "\
                      "(default is coarse-grained)")
    parser.set_defaults(fg=False)
    
    # get the options:
    (options, args) = parser.parse_args()
    
    if not options.pdbfile:
        # if no pdb file was given as option, perhaps as separate arguments:
        if args:
            options.pdbfile = args.pop(0)
        else:
            # if no separate argument remained, we have no input file:
            print("Need at least an input file (pdb)")
            print("")
            parser.print_help()
            print("")
            print("FATAL ERROR: no input file given")
            exit(1)
    if not options.pair:
        # if no chains were given as option, perhaps as separate arguments:
        if args:
            options.chains = args.pop(0)
        else:
            # if no separate argument remained, we have no chains:
            print("Need to specify set of chain IDs for pulling")
            print("")
            parser.print_help()
            print("")
            print("FATAL ERROR: no chain IDs given")
            exit(1)
    if not options.drangefile:
        # if no distances were given as option, perhaps as separate arguments:
        if args:
            options.drangefile = args.pop(0)
        else:
            # if no separate argument remained, we have no distances:
            print("Need to specify distances")
            print("")
            parser.print_help()
            print("")
            print("FATAL ERROR: no distances given")
            exit(1)
            
    # clean up (recommended):
    del(parser)
    return options

#############MAIN##############
if __name__ == "__main__":
    "Do the work"
    options = parse_commandline()

    pdbfile      = options.pdbfile
    pair         = options.pair.split(':')
    drangefile   = options.drangefile
    pull         = options.pull
    boProduction = options.prod
    nsteps       = options.nsteps
    boFG         = options.fg
    
    # check for missing atoms
    if ( pdbt.check_missing_atoms(pdbfile) ):
        # add missing atoms
        mt.add_atoms(pdbfile)
    # minimize original (full atom) pdb structure
    # min.min_orig_pdb(pdbfile)

    # read distances
    try:
        distances = [ d.strip() for line in open(drangefile) ]
    except IOError:
        print("No input file", drangefile, "assuming it is numeric range")
        try:
            start,stop,step = tuple( [ float(dr.strip())
                                       for dr in drangefile.split(":") ] )
            # length of range
            d=stop-start
            # number of steps; including stop (hence the '+1')
            n=int(round(d/step+1.0))
            # calculate nr of decimal places needed for step
            if step>1:
                nd = 0
            else:
                nd = int(abs(math.floor(math.log(0.05, 10))))
            distances = [ str(round(start+i*step, nd)) for i in range(n) ]
        except ValueError:
            print("Not a valid distance range:", drangefile)
            print("Expected 'start:stop:step'; e.g. '2.5:5.0:1.0'")
            exit(1)
    print("Obtained", len(distances), "distances:", distances)

    sim.do_simulations(pdbfile,
                       pair,
                       distances,
                       boProduction,
                       pull,
                       nsteps)

#last line
