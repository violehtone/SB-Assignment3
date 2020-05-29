"""
MFFTools.py

Contains tools transform an atomistic structure to coarse grained structure.

Modified for python3
"""
import os
import re
import sys
import OsTools as ot

import Env as env
import GMX as gmx
import StructTools as st
import PDBTools as pdbt
import GroIO as gio

def ToPdbLineDict(line):
    DictPdbLine            = {}
    DictPdbLine['RecName'] = line[0:7].strip()
    DictPdbLine['Serial']  = line[7:13].strip()
    DictPdbLine['Name']    = line[13:17].strip()
    DictPdbLine['ResName'] = line[17:21].strip()
    DictPdbLine['ChainID'] = line[21].strip()
    DictPdbLine['ResSeq']  = line[22:26].strip()
    DictPdbLine['ICode']   = line[26].strip()
    DictPdbLine['X']       = float(line[30:38].strip())
    DictPdbLine['Y']       = float(line[38:46].strip())
    DictPdbLine['Z']       = float(line[46:55].strip())
    DictPdbLine['Occ']     = float(line[54:60].strip())
    DictPdbLine['TempF']   = float(line[60:66].strip())
    DictPdbLine['Empty']   = ''
    DictPdbLine['Prefix']  = ''
    DictPdbLine['AltLoc']  = ''
#    DictPdbLine['ChainID'] = DictGroups[g]
    DictPdbLine['Element'] = DictPdbLine['Name'][0]
    DictPdbLine['Charge']  = '0'
    return DictPdbLine

def Fg2Cg_NonProtein(PdbFile):
    CgPdbLines = []
    fr = open(PdbFile,'r')
    for line in fr.readlines():
        if line.startswith('ATOM  '):
            PdbLineDict = ToPdbLineDict(line)
            ArgSO4      = (PdbLineDict['ResName']=='SO4' and
                           PdbLineDict['Name']=='S')
            ArgNA       = (PdbLineDict['ResName']=='NA+' and
                           PdbLineDict['Name']=='NA')
            if(ArgSO4 or
               ArgNA):
                CgPdbLines.append(('{0[RecName]:<6}'
                                   '{0[Serial]:>5}'
                                   '{0[Empty]:1}'
                                   '{0[Prefix]:1}'
                                   '{0[Name]:<3}'
                                   '{0[AltLoc]:1}'
                                   '{0[ResName]:<3}'
                                   '{0[Empty]:1}'
                                   '{0[ChainID]:1}'
                                   '{0[ResSeq]:>4}'
                                   '{0[ICode]:1}'
                                   '{0[Empty]:3}'
                                   '{0[X]:>8.3f}'
                                   '{0[Y]:>8.3f}'
                                   '{0[Z]:>8.3f}'
                                   '{0[Occ]:>6.2f}'
                                   '{0[TempF]:>6.2f}'
                                   '{0[Empty]:10}'
                                   '{0[Element]:>2}'
                                   '{0[Charge]:>2}\n'.format(PdbLineDict)))
    fr.close()

    CgPdbFile = ''
    if(len(CgPdbLines)>0):
        CgPdbFile  = re.sub('.pdb','_cg.pdb',PdbFile)
        fw         = open(CgPdbFile,'w')
        for line in CgPdbLines:
            fw.write(line)
        fw.close()

    return CgPdbFile

def finegrained2coarsegrained(fgpair):
    cgpair = []
    for fgfile in fgpair:
        cgfname = re.sub('.pdb','_cg.pdb',fgfile)
        if ot.system(env.fg2cg, options=[fgfile], log=cgfname):
            print("fg2cg failed for", fgfile)
            exit(-1)
        cgpair.append(cgfname)
    return cgpair

def dssp2ssd(dssppair):
    ssdpair = []
    for dsspfile in dssppair:
        ssdname = re.sub('.dssp','.ssd',dsspfile)
        f       = open(ssdname,'w')
        ssd     =  env.to_ssd.dssp2ss(dsspfile)
        f.write(str(len(ssd))+'\n')
        for j in ssd:
            f.write(str(j[1]))
        f.write('\n')
        f.close()
        ssdpair.append(ssdname)
    return ssdpair

def to_itp(fastassddict):
    itppair = []
    for fasta,ssd in fastassddict.items():  # Changed iteritems() to items()
        itpraw = re.sub('.fasta','_raw.itp',fasta)
        itp = re.sub('.fasta','.itp',fasta)
        if ot.system(env.seq2itp,
                     options=['-s '+fasta, '-2 '+ssd, '-t '+itpraw]):
            print("seq2itp failed for", fasta)
            exit(-1)
        fm = open(itpraw,'r').readlines()
        print("Read", len(fm), "from", itpraw)
        molname = fasta.split('.')[0]
        if (len(molname)<8):
            molname += '\t'
        fw = open(itp,'w')
        for line in fm:
            if (re.match('Protein\t1\n',line)):
                line = re.sub('Protein',molname,line)
            fw.write(line)
        fw.close()
        itppair.append(itp)
    return itppair

def select_restr(itpname,ifacelist,elbndsitp):
    ifline = ''
    for i in ifacelist:
        ifline +=' '+i
    if ot.system(env.sel_restr, options=[itpname+ifline], log=elbndsitp):
        print("sel_restr failed for", itpname, ifacelist, elbndsitp)
        exit(-1)

def stabilize_tertiary_structure(groname,
                                 ndxname,
                                 itppair,
                                 groupsuffix,
                                 nat_cgpair):
    stablepair = []
    countitp = 0
    for itp in itppair:
        group     = re.sub('.itp','_'+groupsuffix,itp)
        groupno   = str(st.find_groupnumber(group,ndxname))
        itpname   = group+'_restr.itp'
        if countitp>0:
            # for second group we need postprocessing, write 'raw' file first
            itpraw = group+'_restr_raw.itp'
        else:
            itpraw = itpname
        ifacelist = [ '-f '+groname,
                      '-n '+ndxname,
                      '-constr',
                      '-o '+itpraw
                      ]
        if gmx.g_genrestr(ifacelist, stdin=[groupno],
                          log='genrestr.err'):
            print("Generating restraints failed for", ifacelist)
            exit(-1)
        del ifacelist
#       reindex file if countitp>0
        if (countitp>0):
            itp_mem = open(itpraw,'r').readlines()
            print("Read", len(itp_mem), "from", itpraw)
            f       = open(itpname,'w')
            startindex = 0
            for i in range(len(itp_mem)):
                line = itp_mem[i]
                f.write(line)
                if (re.match('\s*\[\s*constraints\s*\]\s*',line)):
                    startindex = i
                    break
            offset = nat_cgpair[countitp-1]
            for i in range(startindex+1,len(itp_mem)):
                if itp_mem[i].startswith(';'):
                    # silently include comment lines
                    f.write(itp_mem[i])
                else:
                    line = itp_mem[i].split()
                    line[0] = str(int(line[0])-offset)
                    line[1] = str(int(line[1])-offset)
                    line = ' '.join(line)+'\n'
                    f.write(line)
            del itp_mem
            f.close()
        elbndsitp = group+'_addElBnds.itp'
        ifacelist = ['0.5','0.9','500']
        select_restr(itpname, ifacelist, elbndsitp)
        itp_mem = open(itp,'r').readlines()
        print("Read", len(itp_mem), "from", itp)
        lineno  = 0
        for i in range(len(itp_mem)):
            line  = itp_mem[i]
            if (re.match('\s*\[\s*bonds\s*\]\s*',line)):
                lineno = i
                break
        startindex = 0
        for i in range(lineno+1,len(itp_mem)):
            line  = itp_mem[i]
            if (re.match('\[',line)):
                startindex = i-1
                break
        elbndsitp_mem = open(elbndsitp,'r').readlines()
        print("Read", len(elbndsitp_mem), "from", elbndsitp)
        elbndsitp_mem.insert(0,\
          ';elastic bonds needed to stabilize the tertiary structure\n')
        itp_mem[startindex:startindex] = elbndsitp_mem
        del elbndsitp_mem
        fname = re.sub('.itp','_w_elbnds.itp',itp)
        f = open(fname,'w')
        for line in itp_mem:
            f.write(line)
        del itp_mem
        f.close()
        stablepair.append(fname)
        countitp += 1
    return stablepair


def do_fg2cg(pdbfile,pair):

    pdbase  = pdbfile.split('.pdb')[0]

    # file parsing and conversion
    ListChains = []
    if (ot.checkfile(os.getcwd(),pdbfile)):
        splitpdb = 1
        ListChains = env.a2f.run(pdbfile,splitpdb)
    else:
        print('ERROR: checkfile() failed for args \"'+os.getcwd()+\
              '\" and \"'+pdbfile+'\"')
        sys.exit(1)
    
    # First make CG structures of AA chains:
    FastaLs = []
    for file in os.listdir(os.getcwd()):
        if(re.search('fasta',file)):
            FastaLs.append(file)
    DictIsChainAA = {}
    for c in ListChains:
        IsAA = False
        for file in FastaLs:
            if(re.search('_'+c+'.[0-9].fasta',file)):
                IsAA = True
        if(IsAA):
            DictIsChainAA[c] = 'AASequence'
        else:
            DictIsChainAA[c] = 'NonProtein'
    PdbLs = []
    for file in os.listdir(os.getcwd()):
        if(re.search('[0-9].pdb',file)):
            PdbLs.append(file)
            os.system('ls -lfh '+file)
    print("Found pdb files", PdbLs)
    
    FgStructFiles = []
    CgStructFiles = []
    FgFastaFiles  = []
    FgDsspFiles   = []
    FgSsdFiles    = []
    CgItpFiles    = []
    for key,value in DictIsChainAA.items():  # Changed iteritems() to items()
        if(value=='AASequence'):
            FgPdbs   = []
            CgPdbs   = []
            FgFastas = []
            FgDssps  = []
            FgSsds   = []
            CgItps   = []
            for file in PdbLs:
                if(re.search('_'+key+'.[0-9].pdb',file)):
                    FgPdbs.append(file)
                    FgFastas.append(re.sub('.pdb','.fasta',file))
            FgStructFiles.extend(FgPdbs)
            FgFastaFiles.extend(FgFastas)
            for file in FgPdbs:
                CgPdbs.extend(finegrained2coarsegrained([file]))
                FgDssps.extend(pdbt.pdb2dssp([file]))
            for file in FgDssps:
                FgSsds.extend(dssp2ssd([file]))
            DictFastaSsd = {}
            for file in FgSsds:
                FastaFile = re.sub('.ssd','.fasta',file)
                DictFastaSsd[FastaFile] = file
                CgItps.extend(to_itp(DictFastaSsd))
            CgStructFiles.extend(CgPdbs)
            FgDsspFiles.extend(FgDssps)
            FgSsdFiles.extend(FgSsds)
            CgItpFiles.extend(CgItps)
        elif(value=='NonProtein'):
            FgPdbs   = []
            CgPdbs   = []
            CgItps   = []
            for file in PdbLs:
                if(re.search('_'+key+'.[0-9].pdb',file)):
                    FgPdbs.append(file)
                    CgPdbs.append(Fg2Cg_NonProtein(file))
     # print(FgStructFiles
     # print(CgStructFiles
     # print(FgFastaFiles
     # print(FgDsspFiles
     # print(FgSsdFiles
     # print(CgItpFiles

    fgpair    = pdbt.merge_pdbs(pdbase,pair)
    fastapair = pdbt.merge_fastas(pdbase,pair)
    for i in fgpair:
        if (not ot.checkfile(os.getcwd(),i)):
            print('ERROR: checkfile() failed for args \"'+os.getcwd()+\
              '\" and \"'+i+'\"')
            sys.exit(1)
    cgpair   = finegrained2coarsegrained(fgpair)
    dssppair = pdbt.pdb2dssp(fgpair)
    for i in dssppair:
        if (not ot.checkfile(os.getcwd(),i)):
            print('ERROR: checkfile() failed for args \"'+os.getcwd()+\
              '\" and \"'+i+'\"')
            sys.exit(1)
    ssdpair  = dssp2ssd(dssppair)
    itppair  = []
    if (len(ssdpair)==len(fastapair)):
        fastassddict = {}
        i = 0
        for fasta in fastapair:
            fastassddict[fasta] = ssdpair[i]
            i += 1
        for fasta in fastassddict.keys():  # Changed iterkeys() to keys()
            if (not ot.checkfile(os.getcwd(),fasta)):
                print('ERROR: checkfile() failed for args \"'+os.getcwd()+\
                  '\" and \"'+fasta+'\"')
                sys.exit(1)
        for ssd in fastassddict.values():  # Changed itervalues() to values()
            if (not ot.checkfile(os.getcwd(),ssd)):
                print('ERROR: checkfile() failed for args \"'+os.getcwd()+\
                  '\" and \"'+ssd+'\"')
                sys.exit(1)
        itppair = to_itp(fastassddict)
        # check if order of fastapair == order of itppair
        # since dictionary order can be different from list order
        for i in range(len(itppair)):
            itpname = itppair[i]
            itpbase = itpname.split('.')[0]
            itpext  = itpname.split('.')[1]
            fastaname = fastapair[i]
            fastabase = fastaname.split('.')[0]
            if (not re.match(fastabase,itpbase)):
                itppair[i] = fastabase+'.'+itpext
    else:
        print('ERROR: len(ssdpair)!=len(fastapair)')
        sys.exit(1)

    gro_cgpair = gio.to_gro(cgpair)
    nat_cgpair = st.n_atoms(gro_cgpair)
    nrs_cgpair = st.n_residues(gro_cgpair)

    cggro_merged = pdbase+'_cg.gro'
    gio.merge_gro(cggro_merged,gro_cgpair,nat_cgpair,nrs_cgpair)

    # make indexfile
    ifacelist = [ '-f '+cggro_merged,
                  '-o '+'tmp.ndx'
                  ]
    gmx.g_make_ndx(ifacelist, stdin=['q'], log='make_ndx.err')
    del ifacelist
    ngroups = st.ngroups_in_ndx('tmp.ndx')
    os.remove('tmp.ndx')
    cgndx_merged = re.sub('.gro','.ndx',cggro_merged)
    ifacelist = [ '-f '+cggro_merged,
                  '-o '+cgndx_merged
                  ]
    gmx.g_make_ndx(ifacelist,
                   stdin=gmx.gen_ndxUIlist(ngroups,gro_cgpair,nat_cgpair),
                   log='make_ndx.err')
    del ifacelist

    # stabilize the tertiary structure:
    cg_stable_itppair = stabilize_tertiary_structure(cggro_merged,cgndx_merged,itppair,'backbone',nat_cgpair)

    return cggro_merged, cgndx_merged, cg_stable_itppair, \
           gro_cgpair, nat_cgpair, fgpair
