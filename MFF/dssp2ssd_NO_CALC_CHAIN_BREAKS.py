#!/usr/bin/python3
# -*- coding: ISO-8859-1 -*-
#---------------------------------------------------------------------------#
# Function: Extract secondary structure data from a DSSP output file.       #
# Usage: dssp2ssd.py [options] <dssp-output>                                #
# Help: run dssp2ssd.py --help                                              #
# Author: Martti Louhivuori (m.j.louhivuori@rug.nl)                         #
# Version: 0.5 (05.01.2009)                                                 #
#---------------------------------------------------------------------------#
from optparse import OptionParser
import logging, sys, re

def dssp2ss(filename):
    dssp = open(filename)
    ok = False
    data = []
    for line in dssp.readlines():
        if ok and len(line) > 50:
            id = line[6:10]
            ss = line[16]
            aa = line[13:15]
            if ss == ' ':
                ss = '~'
            if aa != '!*':
                if aa != '! ':
                    data.append((id, ss))
        if re.search("#  RESIDUE AA STRUCTURE BP1 BP2", line):
            ok = True
    return data

def echo_ssd(data):
    print(len(data))
    ssd = ''
    for item in data:
        ssd += str(item[1])
    print(ssd)

if __name__ == '__main__':
    usage = 'usage: %prog [options] <dssp-output>'
    desc = 'Extract secondary structure data from a DSSP output file.'
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option('-o', '--output', metavar='FILE', default=None,
            help='write output to FILE (default: STDOUT)')
    parser.add_option('--verbose', action='store_true', default=False,
            help='display additional information while running')
    parser.add_option('--debug', action='store_true', default=False, 
            help='run in debug mode, i.e. maximum information')

    options, args = parser.parse_args()

    # set logger format etc.
    logging.basicConfig(level=logging.WARNING, format='%(levelname)s ' + \
            '%(message)s @ %(asctime)s %(module)s line %(lineno)s',
            datefmt='%H:%M:%S')
    # set logging thresholds
    if options.debug:
        logging.getLogger('').setLevel(logging.DEBUG)
    elif options.verbose:
        logging.getLogger('').setLevel(logging.WARNING)
    else:
        logging.getLogger('').setLevel(logging.CRITICAL)
    logging.debug('options: %s' % repr(options))
    logging.debug('args: %s' % repr(args))

    # redirect STDOUT?
    if options.output:
        try:
            sys.stdout = open(options.output, 'w')
        except IOError:
            print('#', IOError.filename, " - ", IOError.errno)
            print('# Directing output to STDOUT instead.')

    # too few arguments?
    if len(args) < 1:
        parser.error('too few arguments')

    ss = dssp2ss(args[0])
    echo_ssd(ss)

    # the end.
    if options.output:
        sys.stdout.close()

