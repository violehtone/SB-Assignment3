#!/usr/bin/env python3
"""
COMTools.py

Contains tools to change the centre-of-mass distance of two groups of atoms.

Modified for python 3
"""
from math import sqrt
import re
import sys

import GroIO as gio
import OsTools as ot
import StructTools as st

def COM(gro_mem,BBlist):
    atomlist = gio.strip_gro(gro_mem)
    del gro_mem
    r = [0.0,0.0,0.0]
    for i in BBlist:
        llist = atomlist[int(i)-1]
        for j in range(len(r)):
            r[j] += float(llist[j+4])
    del atomlist
    mass = float(len(BBlist))
    for j in range(len(r)):
        r[j] /= max(mass,0.5)
    return r

def distance(ri,rj):
    rij      = []
    distance = 0.0
    for k in range(len(ri)):
        rij.append(rj[k]-ri[k])
        distance += rij[k]**2
    distance = sqrt(distance)
    rij.append(distance)
    return rij

def distance2(ri,rj):
    rij       = []
    distance2 = 0.0
    for k in range(len(ri)):
        rij.append(rj[k]-ri[k])
        distance2 += rij[k]**2
    rij.append(distance2)
    return rij

def com_dist(cggro_merged,cgndx_merged):
    ndx_mem = ot.put_ifile_in_memory(cgndx_merged)
#   make backbone_grouplist
    bbgrouplist = []
    for i in ndx_mem:
        for j in i:
            if (re.search('backbone',j)):
                bbgrouplist.append(j)
    bbatomlist = []
    for group in bbgrouplist:
        bbatomlist.append(st.indexlist(ndx_mem,group))
    del ndx_mem
    gro_mem = gio.put_grofile_in_memory(cggro_merged)
    comlist = []
    for atlist in bbatomlist:
        comlist.append(COM(gro_mem,atlist))
    return distance(comlist[0],comlist[1])

def unit_vector_scalar_product(d,u):
    vec = []
    for i in range(len(u)):
        vec.append(u[i]*d)
    del u
    return vec

def unit_vector(r):
    u = r[0:3]
    for i in range(len(u)):
        u[i] /= r[3]
    del r
    return u

def displace_struct(dcom,groups,d,grofile,ndxfile):
    ndx_mem = ot.put_ifile_in_memory(ndxfile)
    grouplist1 = st.indexlist(ndx_mem,groups[0])
    grouplist2 = st.indexlist(ndx_mem,groups[1])
    del ndx_mem
    del groups
    gro_mem  = gio.put_grofile_in_memory(grofile)
    head_mem = gro_mem[0:2]
    tail_mem = gro_mem[-1]
    coordsnew_mem = gio.strip_gro(gro_mem)
    del gro_mem
    u = unit_vector(dcom)
    d_add = dcom[3]-float(d)
    add_to_gr1 = unit_vector_scalar_product(0.5*d_add,u)
    add_to_gr2 = unit_vector_scalar_product(-0.5*d_add,u)
    for i in grouplist1:
        line = coordsnew_mem[int(i)-1]
        for j in range(len(add_to_gr1)):
            line[4+j] = str(float(line[4+j]) + add_to_gr1[j])
    del grouplist1
    del add_to_gr1
    for i in grouplist2:
        line = coordsnew_mem[int(i)-1]
        for j in range(len(add_to_gr2)):
            line[4+j] = str(float(line[4+j]) + add_to_gr2[j])
    del grouplist2
    del add_to_gr2
    coordsnew_mem[:0] = head_mem
    coordsnew_mem.append(tail_mem)
    del head_mem
    del tail_mem
    ogrofile = re.sub('.gro','_d'+d+'.gro',grofile)
    gio.write_gro(ogrofile,coordsnew_mem)
    return ogrofile

def displace_struct_without_index_file(dcom,grouplist1,grouplist2,d,grofile):
    gro_mem  = gio.put_grofile_in_memory(grofile)
    head_mem = gro_mem[0:2]
    tail_mem = gro_mem[-1]
    coordsnew_mem = gio.strip_gro(gro_mem)
    del gro_mem
    u = unit_vector(dcom)
    d_add = dcom[3]-float(d)
    add_to_gr1 = unit_vector_scalar_product(0.5*d_add,u)
    add_to_gr2 = unit_vector_scalar_product(-0.5*d_add,u)
    for i in grouplist1:
        line = coordsnew_mem[int(i)-1]
        for j in range(len(add_to_gr1)):
            line[4+j] = str(float(line[4+j]) + add_to_gr1[j])
    del grouplist1
    del add_to_gr1
    for i in grouplist2:
        line = coordsnew_mem[int(i)-1]
        for j in range(len(add_to_gr2)):
            line[4+j] = str(float(line[4+j]) + add_to_gr2[j])
    del grouplist2
    del add_to_gr2
    coordsnew_mem[:0] = head_mem
    coordsnew_mem.append(tail_mem)
    del head_mem
    del tail_mem
    ogrofile = re.sub('.gro','_d'+d+'.gro',grofile)
    gio.write_gro(ogrofile,coordsnew_mem)
    return ogrofile

def GenerateBackBoneList(GroMem,istart,iend):
    BBList = []
    for i in range(2+istart-1,2+iend):
        if(re.match('CA',GroMem[i][2]) or
           re.match('BN0',GroMem[i][2])):
            BBList.append(GroMem[i][3])
    return BBList

if __name__ == "__main__":
    "Do the work"
    grofile     = sys.argv[1]
    fr          = open(sys.argv[2],'r')
    FMem        = fr.readlines()
    fr.close()
    istart1     = int(FMem[0].strip().split()[0])
    iend1       = int(FMem[0].strip().split()[1])
#    isolvstart1 = int(FMem[1].strip().split()[0])
#    isolvend1   = int(FMem[1].strip().split()[1])
    del FMem
    fr          = open(sys.argv[3],'r')
    FMem        = fr.readlines()
    fr.close()
    istart2     = int(FMem[0].strip().split()[0])
    iend2       = int(FMem[0].strip().split()[1])
#    isolvstart2 = int(FMem[1].strip().split()[0])
#    isolvend2   = int(FMem[1].strip().split()[1])
    del FMem

    d           = sys.argv[4]

    GroMem        = gio.put_grofile_in_memory(grofile)
    BackBoneList1 = GenerateBackBoneList(GroMem,istart1,iend1)
    BackBoneList2 = GenerateBackBoneList(GroMem,istart2,iend2)
    COM1          = COM(GroMem,BackBoneList1)
    COM2          = COM(GroMem,BackBoneList2)
    DCOM          = distance(COM1,COM2)
    print(DCOM)

    grouplist1 = range(istart1,iend1+1)
#    grouplist1.extend(range(isolvstart1,isolvend1+1))
    grouplist2 = range(istart2,iend2+1)
#    grouplist2.extend(range(isolvstart2,isolvend2+1))
    if(d=='nodisp'):
        d=str(DCOM[3])
    print('Exporting displaced coordinates to file \"'+\
          displace_struct_without_index_file(DCOM,grouplist1,grouplist2,d,grofile)+\
          '\" ...')

