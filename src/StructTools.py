"""
StructTools.py

Contains tools to analyse & modify the structure. E.g. add ions/charges.

Modified for python 3
"""
import sys
import os
import re

import OsTools as ot
import GroIO as gio
import GMX as gmx
import Env as env


def n_residues(gro_in):
    if (type(gro_in)==str):
        if (not ot.checkfile(os.getcwd(),gro_in)):
            print('ERROR (n_residues): checkfile() failed for args \"'+os.getcwd()+\
                  '\" and \"'+gro_in+'\"')
            sys.exit(1)
        gro_mem = gio.put_grofile_in_memory(gro_in)
        nrs = int(gro_mem[-2][0])
        del gro_mem
        return nrs
    elif (type(gro_in)==list):
        nrs_list = []
        for file in gro_in:
            if (not ot.checkfile(os.getcwd(),file)):
                print('ERROR (n_residues): checkfile() failed for args \"'+os.getcwd()+\
                      '\" and \"'+file+'\"')
                sys.exit(1)
            gro_mem = gio.put_grofile_in_memory(file)
            nrs     = int(gro_mem[-2][0])
            del gro_mem
            nrs_list.append(nrs)
        return nrs_list


def n_atoms(gro_in):
    if (type(gro_in)==str):
        if (not ot.checkfile(os.getcwd(),gro_in)):
            print('ERROR (n_atoms): checkfile() failed for args \"'+os.getcwd()+\
                  '\" and \"'+gro_in+'\"')
            sys.exit(1)
        gro_mem = gio.put_grofile_in_memory(gro_in)
        nat = int(gro_mem[1][0])
        del gro_mem
        return nat
    elif (type(gro_in)==list):
        nat_list = []
        for file in gro_in:
            if (not ot.checkfile(os.getcwd(),file)):
                print('ERROR (n_atoms): checkfile() failed for args \"'+os.getcwd()+\
                      '\" and \"'+file+'\"')
                sys.exit(1)
            gro_mem = gio.put_grofile_in_memory(file)
            nat     = int(gro_mem[1][0])
            del gro_mem
            nat_list.append(nat)
        return nat_list


def ngroups_in_ndx(ndx_in):
    f = open(ndx_in,'r')
    content = f.read()
    m = re.findall('\[',content)
    f.close()
    del content
    return len(m)


def indexlist(fi_mem,m_expr):
    rr        = re.compile("\[")
    readstart = 0
    for i in range(len(fi_mem)):
        line = fi_mem[i]
        for j in range(len(line)):
            if(rr.match(line[j])):
                group = line[j+1]
                if(group==m_expr):
                    readstart = i
    grouplist = []
    for i in range(readstart+1,len(fi_mem)):
        linelist = fi_mem[i]
        if(len(linelist)>0):
            if(linelist[0]=='['):
                break
        for j in linelist:
            grouplist.append(j)
    return grouplist


def mindist_between_groups(cggrofile,cgndxfile,groups):
    mindistfile = 'mindist.xvg'
    ifacelist = [ '-f '+cggrofile,
                  '-n '+cgndxfile
                  ]
    stdin = [ str(find_groupnumber(groups[0],cgndxfile)),
              str(find_groupnumber(groups[1],cgndxfile))
              ]
    gmx.g_g_mindist(ifacelist, stdin=stdin, log='mindist.err ')
    mindist_mem = ot.put_ifile_in_memory(mindistfile)
    mindist = float(mindist_mem[-1][-1])
    del mindist_mem
    return mindist
#    ndx_mem = put_ifile_in_memory(cgndx_merged)
##   make backbone_grouplist
#    atomlist = []
#    for group in groups:
#        atomlist.append(indexlist(ndx_mem,group))
#    del ndx_mem
#    gro_mem = put_grofile_in_memory(cggro_merged)
#    at_mem  = strip_gro(gro_mem)
#    del gro_mem
#    ri = [0.0,0.0,0.0]
#    rj = [0.0,0.0,0.0]
#    mindist = 1.0e100
#    for listA in range(len(atomlist)-1):
#        for listB in range(listA+1,len(atomlist)):
#            for atomAi in atomlist[listA]:
#                lA = at_mem[int(atomAi)-1]
#                for i in range(len(ri)):
#                    ri[i] = float(lA[i+4])
#                for atomBj in atomlist[listB]:
#                    lB = at_mem[int(atomBj)-1]
#                    for j in range(len(rj)):
#                        rj[j] = float(lB[j+4])
#                    mindist = min(mindist,distance2(ri,rj)[3])
#    mindist = sqrt(mindist)
#    return mindist


def find_groupnumber(group,ndxname):
    ndx_mem = ot.put_ifile_in_memory(ndxname)
    countgroup = 0
    for line in ndx_mem:
        if ((len(line)>0) and (re.match('\[',line[0]))):
            if (re.match(group,line[1])):
                break
            countgroup += 1
    del ndx_mem
    return countgroup


def get_charge_from_itp(itp,nat):
    f = open(itp,'r')
    itp_mem = f.readlines()
    f.close()
    startindex = 0
    for i in range(len(itp_mem)):
        line = itp_mem[i]
        if(re.search('\[atoms\]',line)):
            startindex = i+1
            break
    q = 0
    for i in range(nat):
        index = i + startindex
        l = itp_mem[index].strip().split()
        try:
            q += int(float(l[6]))
        except ValueError:
            print("cannot parse value", l[6], "from line", index, "in file", itp)
            exit(-1)
    del itp_mem
    return q


def neutralize(cggro,itppair,nat_cgpair):
    q = 0
    for i in range(len(itppair)):
        itp = itppair[i]
        nat = nat_cgpair[i]
        q  += get_charge_from_itp(itp,nat)
    cggro_neutral = re.sub('.gro','_neutral.gro',cggro)
    if (q==0):
#       add nothing
        if ot.cp(cggro, cggro_neutral):
            print("cp failed for", cggro, cggro_neutral)
            exit(-1)
    else:
        ifacelist = []
        ifacelist.append('-cp '+cggro)
        if (q>0):
#           add chloride
            ifacelist.append('-ci '+env.mffstruct.replace(' ','\ ')+'/cl-.gro')
            ifacelist.append('-nmol '+str(q))
        elif (q<0):
#           add sodium
            q = abs(q)
            ifacelist.append('-ci '+env.mffstruct.replace(' ','\ ')+'/na+.gro')
            ifacelist.append('-nmol '+str(q))
        ifacelist.append('-try 1000')
        ifacelist.append('-vdwd 0.19')
        ifacelist.append('-o '+cggro_neutral)
        gmx.g_genbox(ifacelist, log='genbox.err')
        del ifacelist
    return cggro_neutral


def gen_posre(ndxfile,grofile,groups,nat_cgpair):
    posrelist = []
    for i in range(len(groups)):
        ifacelist=['-f '+grofile]
        ifacelist.append('-n '+ndxfile)
        posreitp = groups[i]+'_posre.itp'
        if (i>0):
            # for second group we need postprocessing, write 'raw' file first
            posreraw = groups[i]+'_posre_raw.itp'
        else:
            posreraw = posreitp
        posrelist.append(posreitp)
        ifacelist.append('-o '+posreraw)
        groupno = str(find_groupnumber(groups[i],ndxfile))
        gmx.g_genrestr(ifacelist, stdin=[groupno], log='genrestr.err')
        if (i>0):
            posre_mem = open(posreraw,'r').readlines()
            fw = open(posreitp,'w')
            startindex = 0
            for j in range(len(posre_mem)):
                line = posre_mem[j]
                fw.write(line)
                if (re.match(';  i funct       fcx        fcy        fcz',line)):
                    startindex = j
                    break
            offset = nat_cgpair[i-1]
            for j in range(startindex+1,len(posre_mem)):
                line    = posre_mem[j].split()
                line[0] = str(int(line[0])-offset)
                fw.write('%4s%6s%10s%11s%11s\n' % (line[0],line[1],line[2],line[3],line[4]))
            fw.close()
            del posre_mem
    return posrelist
