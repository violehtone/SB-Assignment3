"""
GroIO.py

Contains tools for handling gromacs coordinate files (.gro).

Modified for python 3
"""
import sys
import os
import re

import OsTools as ot
import GMX as gmx

def put_grofile_in_memory(ifile):
    fi     = open(ifile,"r")
    fi_mem = fi.readlines()
    fi.close()
    length = len(fi_mem)
    for i in range(2):
        l         = fi_mem[i].strip().split()
        fi_mem[i] = l
    for i in range(2,length-1):
        l  = fi_mem[i]
        ll = []
#       .gro format:
#       First 4 fields of five columns of .gro format
#       Residue number, Residue type, Atom type, Atom number
        count = 0
        for field in range(4):
            start = field*5
            end   = start+5
            count = end
            ll.append(l[start:end].strip())
#       Fields 4 to 7 of 8 columns of .gro format
#       x-, y-, z-coordinate
        startcount = count
        for field in range(3):
            start = startcount + field*8
            end   = start + 8
            count = end
            ll.append(l[start:end].strip())
#       Fields 7 to 10 of .gro format
#       x-, y-, z-velocities
        if (count < len(l)-2):
            start = count
            end   = len(l)-1
            ll.extend(l[start:end].strip().split())
        fi_mem[i] = ll
        del ll
    for i in range(length-1,length):
        l         = fi_mem[i].strip().split()
        fi_mem[i] = l

    return fi_mem

def to_gro(pdb_in):
    if (type(pdb_in)==str):
        if (not ot.checkfile(os.getcwd(),pdb_in)):
            print('ERROR (to_gro): checkfile() failed for args \"'+os.getcwd()+\
                  '\" and \"'+pdb_in+'\"')
            sys.exit(1)
        gro_out = re.sub('.pdb','.gro',pdb_in)
        ifacelist = [ '-f '+pdb_in,
                      '-o '+gro_out
                      ]
        gmx.g_editconf(ifacelist, log='editconf.err')
        del ifacelist
        return gro_out
    elif (type(pdb_in)==list):
        gro_out_list = []
        for file in pdb_in:
            if (not ot.checkfile(os.getcwd(),file)):
                print('ERROR (to_gro): checkfile() failed for args \"'+os.getcwd()+\
                      '\" and \"'+file+'\"')
                sys.exit(1)
            gro_out = re.sub('.pdb','.gro',file)
            ifacelist = [ '-f '+file,
                          '-o '+gro_out
                          ]
            gmx.g_editconf(ifacelist, log='editconf.err')
            del ifacelist
            gro_out_list.append(gro_out)
        return gro_out_list

def strip_gro(gro_mem):
    gro_stripped = gro_mem[2:-1]
    del gro_mem
    return gro_stripped

def write_gro(ofile,new_coords_mem):
    f = open(ofile,"w")
    new_coords_mem.reverse()
    f.write("%s\n" % ''.join(new_coords_mem.pop()))
    natom = ''.join(new_coords_mem.pop())
    f.write("%s\n" % natom)
    natom = int(natom)
    for i in range(0,natom):
        l = new_coords_mem.pop()
        if (i<9999):
            rx = float(l[4])
            ry = float(l[5])
            rz = float(l[6])
            if(len(l)>7):
                vx = float(l[7])
                vy = float(l[8])
                vz = float(l[9])
                f.write("%5s%5s%5s%5s%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n" % \
                (l[0],l[1],l[2],l[3],rx,ry,rz,vx, vy, vz))
            else:
                f.write("%5s%5s%5s%5s%8.3f%8.3f%8.3f\n" % \
                (l[0],l[1],l[2],l[3],rx,ry,rz))
#        else:
#            rx = float(lsplit[2])
#            ry = float(lsplit[3])
#            rz = float(lsplit[4])
#            if(len(lsplit)>5):
#                vx = float(lsplit[5])
#                vy = float(lsplit[6])
#                vz = float(lsplit[7])
#                f.write("%8s %11s %7.3f %7.3f %7.3f %7.4f %7.4f %7.4f\n" % \
#                (lsplit[0],lsplit[1],rx,ry,rz,vx, vy, vz))
#            else:
#                f.write("%8s %11s %7.3f %7.3f %7.3f\n" % \
#                (lsplit[0],lsplit[1],rx,ry,rz))
    f.write("  %s\n" % '   '.join(new_coords_mem.pop()))
    f.close()
    del new_coords_mem

def merge_gro(fname,gro_cgpair,nat_cgpair,nrs_cgpair):
#   use only with protein residues
    gro_merged_mem = []
    title = fname.split('.')[0]
    gro_merged_mem.append([title])
    nat_tot = 0
    for i in nat_cgpair:
        nat_tot += i
    gro_merged_mem.append([str(nat_tot)])
    box = ""
    for file in range(len(gro_cgpair)):
        tmp_mem = put_grofile_in_memory(gro_cgpair[file])
        gro_mem = strip_gro(tmp_mem)
        box     = tmp_mem[-1]
        del tmp_mem
        res_offset = 0
        at_offset  = 0
        if (file>0):
            res_offset = nrs_cgpair[file-1]
            at_offset  = nat_cgpair[file-1]
        for line in gro_mem:
            line[0] = str(int(line[0])+res_offset)
            line[3] = str(int(line[3])+at_offset)
            gro_merged_mem.append(line)
        del gro_mem
    gro_merged_mem.append(box)
    write_gro(fname,gro_merged_mem)
    del gro_merged_mem
