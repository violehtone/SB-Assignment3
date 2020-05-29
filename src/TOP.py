"""
TOP.py

Contains tools for generating MARTINI topology files (.top).

Modified for python 3
"""
import os
import re

import GroIO as gio

def ChangeNPartInTop(N,MolType):
    FileBase = 'cg_'+MolType+str(N)
    return generate_top_general(filebase=FileBase,Type=MolType,NPart=N)


def generate_top_general(filebase='cg',sysname='water',grofile='EMPTY',Type='',NPart=0):
    ext = 'top'
    topfile = filebase+'.'+ext
    f = open(topfile,'w')
    f.write(';\n')
    f.write(';'+filebase+' | MARTINI 2.1\n')
    f.write(';\n')
    f.write('\n')
    f.write('; Force field:\n')
    ls = os.listdir('./')
    for dirfile in ls:
        if(re.search('martini_v2.1.itp',dirfile)):
            f.write('#include \"'+dirfile+'\"\n')
            ls.remove(dirfile)
    for dirfile in ls:
        if(re.search('martini_v2.0_ions.itp',dirfile)):
            f.write('#include \"'+dirfile+'\"\n')
            ls.remove(dirfile)
    f.write('\n')
    f.write('[ system ]\n')
    f.write('; Name\n')
    f.write(sysname+'\n')
    f.write('\n')
    f.write('[ molecules ]\n')
    f.write('; Compound\tNmol\n')
    if(grofile!='EMPTY'):
        gro_mem = gio.put_grofile_in_memory(grofile)
        for solvent_molecule in ['NA+','CL-','W']:
            count = 0
            for line in gro_mem:
                if ((len(line)>2) and (re.match(solvent_molecule,line[2]))):
                    count += 1
        if (count>0):
            name = solvent_molecule
        if (len(name)<8):
            name += '\t'
        f.write(name+'\t'+str(count)+'\n')
        del gro_mem
    else:
        if (len(Type)<8):
            Type += '\t'
        f.write(Type+'\t'+str(NPart)+'\n')
    f.close
    return topfile


def generate_top(filebase,pdbid,itppair,fgpair,grofile,posrepair=[]):
    ext = 'top'
    topfile = filebase+'.'+ext
    f = open(topfile,'w')
    f.write(';\n')
    f.write(';'+filebase+' | MARTINI 2.1\n')
    f.write(';\n')
    f.write('\n')
    f.write('; Force field:\n')
    ls = os.listdir('./')
    print('Including martini*.itp')
    for dirfile in ls:
        if(re.search('martini_v2.1.itp',dirfile)):
            f.write('#include \"'+dirfile+'\"\n')
            ls.remove(dirfile)
    for dirfile in ls:
        if(re.search('martini_v2.0_ions.itp',dirfile)):
            f.write('#include \"'+dirfile+'\"\n')
            ls.remove(dirfile)
    f.write('\n')
    if (filebase=='cg' or filebase=='cg_posre'):
        print('Including', itppair, posrepair)
        if not posrepair:
            posrepair = [ None for i in itppair ]
        for itpfile,posrefile in zip(itppair,posrepair):
            for dirfile in sorted(ls):
                if(re.search(itpfile,dirfile)):
                    f.write('#include \"'+dirfile+'\"\n')
                    ls.remove(dirfile)
            for dirfile in sorted(ls):
                if posrefile:
                    if(re.search(posrefile,dirfile)):
                        f.write('#ifdef POSRES\n')
                        f.write('#include \"'+posrefile+'\"\n')
                        f.write('#endif\n')
                        ls.remove(dirfile)
    f.write('\n')
    f.write('[ system ]\n')
    f.write('; Name\n')
    f.write(pdbid+'\n')
    f.write('\n')
    f.write('[ molecules ]\n')
    f.write('; Compound\tNmol\n')
    for name in fgpair:
        molname = name.split('.')[0]
        if (len(molname)<8):
            molname += '\t'
        f.write(molname+'\t1\n')
    gro_mem = gio.put_grofile_in_memory(grofile)
    for solvent_molecule in ['NA+','CL-','W']:
        count = 0
        for line in gro_mem:
            if ((len(line)>2) and (re.match(solvent_molecule,line[2]))):
                count += 1
        if (count>0):
            name = solvent_molecule
            if (len(name)<8):
                name += '\t'
            f.write(name+'\t'+str(count)+'\n')
    del gro_mem
    f.close
    return topfile
