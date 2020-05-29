#!/usr/bin/python3
"""
Modified for python3
"""
import os
import re
import sys

import OsTools as ot
import COMTools as ct
import StructTools as st
import GMX as gmx
import MDP as mdp
import Env as env
import TOP as top

def do_one_sim(name='em_vac',  # name of run part; used as dirname
               run_type='em',  # type of run: em, em_posre, 
               		       # md_pull_eq, md_pull_prod, md_umbr_prod
               top_type='cg',  # type of topology: cg, cg_posre
               prev_conf='../conf.gro', # input conformation
               prev_ndx =None, # input index file
               pdbase='system',# name for system/molecule in topology file
               cg_stable_itppair='', # some itp file?
               fgpair='',      # some itp file?
               posrepair=[],   # position restraints itp file,
               posredir=None,  # dir where to find pr itp file
               groups='',      # index groups for chains across interface
               nat_cgpair=[],  # nr of atoms in each chain
               sdistance=None, # con/restraint distance
               boRun=True      # start mdrun
               ):
    # create and change to run dir:
    prev_dir = os.getcwd()
    dir = name
    ot.GotoDir(dir)
    
    # link input conformation file:
    this_conf = 'conf.gro'
    ot.symlink(prev_conf,this_conf)

    # generate MD parameter file:
    if sdistance:
        # add pull groups to MDP:
        bbgroups  = [ group+'_backbone' for group in groups ]
    else:
        bbgroups = []
    mdpfile = mdp.generate_mdp(run_type, sdistance, bbgroups)
    
    # we rely on being two dir levels deep: e.g., 'd2.30/md_sol_pr'
    for file in cg_stable_itppair:
        ot.symlink('../../'+file,'./'+file)
    
    # link to Martini base topology files:
    if ot.symlink(env.mffFF+'/'+env.mffitp,  env.mffitp):  exit(-1)
    if ot.symlink(env.mffFF+'/'+env.mffions, env.mffions): exit(-1)
    
    if prev_ndx:
        # link index file
        this_ndx=os.path.split(prev_ndx)[1] # returns ['/path', 'file']
        ot.symlink(prev_ndx,this_ndx)
    else:
        this_ndx=None
    
    ret_posrepair = []
    if posrepair:
        # make local links:
        this_posrepair=[]
        for file in posrepair:
            this = os.path.split(file)[1] # returns ['/path', 'file']
            ot.symlink('../'+posredir+'/'+file, this)
            this_posrepair.append(this)
        posrepair=this_posrepair
        # now trim pathname from posrepair file names:
    elif this_ndx:
        # generate posre files:
        print("Gen POSRE", this_ndx, this_conf, groups, nat_cgpair)
        posrepair = st.gen_posre(this_ndx, this_conf, groups, nat_cgpair)
        print("Got POSRE", posrepair)
        ret_posrepair = posrepair
    topfile = top.generate_top(top_type, pdbase, cg_stable_itppair,
			       fgpair, this_conf, posrepair)
    
    # grompp
    optionlist = [ '-f '+mdpfile,
                   '-c '+this_conf,
                   '-p '+topfile
                   ]
    if this_ndx:                    optionlist.append('-n '+this_ndx)
    if run_type.startswith('em'):   optionlist.append('-maxwarn 1')
    elif run_type.startswith('md'): optionlist.append('-maxwarn 3')
    result = gmx.g_grompp(optionlist, log='grompp.err')
    # return immediately if error occurred
    if result: return result, dir, ret_posrepair
    
    # mdrun
    if boRun:
        optionlist = [ '-s topol.tpr' ]
        result = gmx.g_mdrun(optionlist, log='md.err')
    else:
        print('Ready for production in', os.getcwd())
        print('option list', optionlist)

    # return to original (parent) dir
    os.chdir(prev_dir)
    
    # next step needs to know our dirname, and posre files
    return result, dir, ret_posrepair

def do_solvate(name='solvate', # name of run part; used as dirname
               prev_conf='../conf.gro',  # input conformation
               prev_ndx ='../index.ndx', # input index file
               gro_cgpair='',  # gro files?
               nat_cgpair=0    # nr of atoms?
               ):
    # solvate structure with cg water
    # create and change to run dir:
    prev_dir = os.getcwd()
    dir = name
    ot.GotoDir(dir)
    
    # link input conformation file:
    this_conf = 'conf.gro'
    ot.symlink(prev_conf,this_conf)

    # solvate structure with cg water
    optionlist = [ '-cp '+this_conf,
                   '-cs '+env.mffstruct.replace(' ','\ ')+'/water.gro',
                   '-vdwd 0.19'
                   ]
    result = gmx.g_genbox(optionlist, log='genbox.err')
    if result: return result
    
    # make a new ndx file of the solvated structure
    outgro = 'out.gro'
    tmpndx = 'tmp.ndx'
    optionlist = [ '-f '+outgro,
                   '-o '+tmpndx
                   ]
    stdin =      [ 'q' ] # just quit; writes out all default index groups
    result = gmx.g_make_ndx(optionlist, log='make_ndx.err', stdin=stdin)
    if result: return result
    ngroups = st.ngroups_in_ndx(tmpndx)
    os.remove(tmpndx)
    
    ndx_sol = re.sub('.ndx','_sol.ndx',prev_ndx)
    ndxlist = gmx.gen_ndxUIlist(ngroups,gro_cgpair,nat_cgpair)
    optionlist = [ '-f '+outgro,
                   '-o '+ndx_sol,
                   ]
    result = gmx.g_make_ndx(optionlist, log='make_ndx.err', stdin=ndxlist)
    if result: return result
    
    # return to original (parent) dir
    os.chdir(prev_dir)
    
    # next step needs to know our dirname and output index and gro files:
    return result, dir, ndx_sol, outgro
    

def do_pullsim(root_dir,distance,cggro_merged,cgndx_merged,
               groups,cg_stable_itppair,
               nat_cgpair,pdbase,fgpair,gro_cgpair,
               boRun=True, pull='cons', nsteps=None):

    os.chdir(root_dir)

    # check if distance is string or float:
    try:
        distance += 0.0
        # is float, make string version
        sdistance = str(distance)
    except TypeError:
        # is not float
        try:
            sdistance = distance+''
            distance = float(sdistance)
        except TypeError:
            print('Unknown type for distance:', distance, type(distance))
                
    # make dir for this distance and cd to dir
    d_dir = os.getcwd()+'/d'+sdistance
    ot.GotoDir(d_dir)
    ot.symlink(root_dir+'/'+cggro_merged,'./'+cggro_merged)
    ot.symlink(root_dir+'/'+cgndx_merged,'./'+cgndx_merged)
    
    # displace structure
    dcom_of_cgpdb = ct.com_dist(cggro_merged,cgndx_merged)
    cggro_disp = ct.displace_struct(dcom_of_cgpdb,groups,sdistance,
                                    cggro_merged,cgndx_merged)

    # center structure in box
    cggro_cntr = re.sub('.gro','_cntr.gro',cggro_disp)
    # distance to box edge should be more than mindist+0.5*rcut (rcut = 1.2 nm)
    distance_to_box_edge = 1.0+st.mindist_between_groups(cggro_merged,
                                                         cgndx_merged,
                                                         groups)
    optionlist = [ '-f '+cggro_disp,
                   '-o '+cggro_cntr,
                   # use default rectangular/triclinic; alternate box shapes:
                   # '-bt cubic', 
                   '-bt dodecahedron',
                   '-c',
                   '-princ',
                   '-d '+str(distance_to_box_edge)
                   ]
    result = gmx.g_editconf(optionlist, log='editconf.err', stdin=['0'])
    if result: return result

    # neutralize the total structure
    for itp in cg_stable_itppair:
        ot.symlink(root_dir+'/'+itp,'./'+itp)
    cggro_neutral = st.neutralize(cggro_cntr,cg_stable_itppair,nat_cgpair)

    # vacuum energy minimization
    result, em_vac_dir, ignore = \
        do_one_sim(name='em_vac',
                   run_type='em', 
                   top_type='cg',
                   prev_conf='../'+cggro_neutral,
                   pdbase=pdbase,
                   cg_stable_itppair=cg_stable_itppair,
                   fgpair=fgpair
                   )
    if result: return result

    # vacuum energy minimization using position restraints
    result, em_vac_pr_dir, posrepair = \
        do_one_sim(name='em_vac_posre',
                   run_type='em_posre', 
                   top_type='cg_posre',
                   prev_conf='../'+em_vac_dir+'/confout.gro',
                   prev_ndx ='../../'+cgndx_merged,
                   pdbase=pdbase,
                   cg_stable_itppair=cg_stable_itppair,
                   fgpair=fgpair,
                   groups=groups,
                   nat_cgpair=nat_cgpair
                   )
    if result: return result

    # solvate structure with cg water
    result, solvate_dir, cgndx_sol, outgro = \
        do_solvate(name='solvate',
                   prev_conf  = '../'+em_vac_pr_dir+'/conf.gro',
                   prev_ndx   = cgndx_merged,
                   gro_cgpair = gro_cgpair,
                   nat_cgpair = nat_cgpair
                   )
    if result: return result

    # solution energy minimization using position restraints
    result, em_sol_pr_dir, ignore = \
        do_one_sim(name='em_sol_posre',
                   run_type='em_posre', 
                   top_type='cg_posre',
                   prev_conf='../'+solvate_dir+'/'+outgro,
                   prev_ndx ='../'+solvate_dir+'/'+cgndx_sol,
                   pdbase=pdbase,
                   cg_stable_itppair=cg_stable_itppair,
                   fgpair=fgpair,
                   posrepair=posrepair,
                   posredir=em_vac_pr_dir
                   )
    if result: return result

    # solution md using position restraints
    result, md_sol_pr_dir, ignore = \
        do_one_sim(name='md_sol_posre',
                   run_type='md_posre', 
                   top_type='cg_posre',
                   prev_conf='../'+em_sol_pr_dir+'/confout.gro',
                   prev_ndx ='../'+solvate_dir+'/'+cgndx_sol,
                   pdbase=pdbase,
                   cg_stable_itppair=cg_stable_itppair,
                   fgpair=fgpair,
                   posrepair=posrepair,
                   posredir=em_vac_pr_dir
                   )
    if result: return result

    # solution equilibration md run
    result, md_sol_eq_dir, ignore = \
        do_one_sim(name='md_sol_eq',
                   run_type='md_pull_eq', 
                   top_type='cg',
                   prev_conf='../'+md_sol_pr_dir+'/confout.gro',
                   prev_ndx ='../'+solvate_dir+'/'+cgndx_sol,
                   pdbase=pdbase,
                   cg_stable_itppair=cg_stable_itppair,
                   fgpair=fgpair,
                   groups=groups,
                   nat_cgpair=nat_cgpair,
                   sdistance=sdistance
                   )
    if result: return result

    # solution production md run
    if   pull=='cons': run='md_pull_prod'
    elif pull=='umbr': run='md_umbr_prod'
    else:
        print("Unknown pull option", pull)
        exit(-1)
    result, md_sol_prod_dir, ignore = \
        do_one_sim(name='md_sol_prod',
                   run_type=run, 
                   top_type='cg',
                   prev_conf='../'+md_sol_eq_dir+'/confout.gro',
                   prev_ndx ='../'+solvate_dir+'/'+cgndx_sol,
                   pdbase=pdbase,
                   cg_stable_itppair=cg_stable_itppair,
                   fgpair=fgpair,
                   groups=groups,
                   nat_cgpair=nat_cgpair,
                   sdistance=sdistance,
                   boRun=boRun
                   )
    if result: return result
    os.chdir(d_dir)

if __name__ == "__main__":
    "Do the work"
    root_dir          = sys.argv[1]
    distance          = sys.argv[2]
    cggro_merged      = sys.argv[3]
    cgndx_merged      = sys.argv[4]
    dcom_of_cgpdb     = [sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8]]
    for i in range(len(dcom_of_cgpdb)):
        fi = float(dcom_of_cgpdb[i])
        dcom_of_cgpdb[i] = fi
    groups            = [sys.argv[9],sys.argv[10]]
    cg_stable_itppair = [sys.argv[11],sys.argv[12]]
    nat_cgpair        = [sys.argv[13],sys.argv[14]]
    for i in range(len(nat_cgpair)):
        ii = int(nat_cgpair[i])
        nat_cgpair[i] = ii
    pdbase            = sys.argv[15]
    fgpair            = [sys.argv[16],sys.argv[17]]
    gro_cgpair        = [sys.argv[18],sys.argv[19]]

    do_pullsim(root_dir,distance,cggro_merged,cgndx_merged,
               dcom_of_cgpdb,groups,cg_stable_itppair,
               nat_cgpair,pdbase,fgpair,gro_cgpair)
