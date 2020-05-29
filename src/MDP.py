"""
MDP.py

Tools to generate MDP files for the Pull simulations.

Modified for python 3
"""
import os
import datetime as dt

def default_preprocess_dict_mdp(type):
    dict = {}
    if (type=='md_pull_eq' or type=='md_pull_prod' or type=='md_umbr_prod'):
        dict = {
                'comment':'; Preprocessing options:',
                'title':type,
                'cpp':'/lib/cpp'
                }
    if (type=='md_posre'):
        dict = {
                'comment':'; Preprocessing options:',
                'title':type,
                'cpp':'/lib/cpp',
                'define':'-DPOSRES'
                }
    if (type=='em_posre'):
        dict = {
                'comment':'; Preprocessing options:',
                'title':type,
                'cpp':'/lib/cpp',
                'define':'-DPOSRES'
                }
    if (type=='em'):
        dict = {
                'comment':'; Preprocessing options:',
                'title':type,
                'cpp':'/lib/cpp'
                }
    if (type=='md'):
        dict = {
                'comment':'; Preprocessing options:',
                'title':type,
                'cpp':'/lib/cpp'
                }
    if (type=='md_single_point'):
        dict = {
                'comment':'; Preprocessing options:',
                'title':type,
                'cpp':'/lib/cpp'
                }
    if (type=='md_move'):
        dict = {
                'comment':'; Preprocessing options:',
                'title':type,
                'cpp':'/lib/cpp'
                }
    if (type=='hybrid'):
        dict = {
                'comment':'; Preprocessing options:',
                'title':type,
                'cpp':'/lib/cpp'
                }
    if (type=='hybrid_sp'):
        dict = {
                'comment':'; Preprocessing options:',
                'title':type,
                'cpp':'/lib/cpp'
                }
    return dict

def default_run_control_dict_mdp(type,nsteps):
    dict = {}
    if (type=='md_single_point'):
        dict = {
                'comment':'; Run control parameters:',
                'integrator':'md',
                'tinit':'0.0',
                'dt':'0.02',
                'nsteps':'0',
                'nstcomm':'1',
                'comm_grps':'System'
                }
    if (type=='md'):
        dict = {
                'comment':'; Run control parameters:',
                'integrator':'md',
                'tinit':'0.0',
                'dt':'0.02',
                'nsteps':str(nsteps),
                'nstcomm':'1',
                'comm_grps':'System'
                }
    if (type=='hybrid'):
        dict = {
                'comment':'; Run control parameters:',
                'integrator':'md',
                'tinit':'0.0',
                'dt':'0.02',
                'nsteps':str(nsteps),
                'nstcomm':'1',
                'comm_grps':'System'
                }
    if (type=='hybrid_sp'):
        dict = {
                'comment':'; Run control parameters:',
                'integrator':'md',
                'tinit':'0.0',
                'dt':'0.02',
                'nsteps':'0',
                'nstcomm':'1',
                'comm_grps':'System'
                }
    if (type=='md_move'):
        dict = {
                'comment':'; Run control parameters:',
                'integrator':'md',
                'tinit':'0.0',
                'dt':'0.02',
                'nsteps':str(nsteps),
                'nstcomm':'1',
                'comm_grps':'System'
                }
    if (type=='md_pull_prod' or type=='md_umbr_prod'):
        dict = {
                'comment':'; Run control parameters:',
                'integrator':'md',
                'tinit':'0.0',
                'dt':'0.02',
                'nsteps':'25000',
                'nstcomm':'1',
                'comm_grps':'System'
                }
    if (type=='md_pull_eq'):
        dict = {
                'comment':'; Run control parameters:',
                'integrator':'md',
                'tinit':'0.0',
                'dt':'0.02',
                'nsteps':'2000',
                'nstcomm':'1',
                'comm_grps':'System'
                }
    if (type=='md_posre'):
        dict = {
                'comment':'; Run control parameters:',
                'integrator':'md',
                'tinit':'0.0',
                'dt':'0.02',
                'nsteps':'2000'
                }
    if (type=='em_posre' or type=='em'):
        dict = {
                'comment':'; Run control parameters:',
                'integrator':'steep',
                'tinit':'0.0',
                'dt':'0.002',
                'nsteps':'10000',
                'em_tol':'550',
                'emstep':'0.01',
                'nstcomm':'1',
                'comm_grps':'System'
                }
    return dict

def default_output_ctrl_dict_mdp(type,nsteps,extragroups):
    dict = {}
    if (type=='em' or type=='em_posre'):
        dict = {
                'comment':'; Output control options:',
                'nstxout':'100',
                'nstvout':'100',
                'nstlog':'100',
                'nstenergy':'100',
                'nstxtcout':'100',
                'energygrps':'System'
                }
    if (type=='md_posre'):
        dict = {
                'comment':'; Output control options:',
                'nstxout':'1000',
                'nstvout':'1000',
                'nstlog':'1000',
                'nstenergy':'1000',
                'nstxtcout':'1000',
                'energygrps':'Protein Non-Protein'
                }
    if (type=='md_pull_eq'):
        dict = {
                'comment':'; Output control options:',
                'nstxout':'1000',
                'nstvout':'1000',
                'nstfout':'1000',
                'nstlog':'1000',
                'nstenergy':'1000',
                'nstxtcout':'1000',
                'xtc_precision':'1000',
                'energygrps':'Protein Non-Protein'
                }
    if (type=='md_pull_prod' or type=='md_umbr_prod'):
        dict = {
                'comment':'; Output control options:',
                'nstxout':'10000',
                'nstvout':'10000',
                'nstfout':'10000',
                'nstlog':'10000',
                'nstenergy':'5',
                'nstxtcout':'10000',
                'xtc_precision':'1000',
                'energygrps':'Protein Non-Protein'
                }
    nrgygrps = ''
    if(len(extragroups)>0):
        for line in extragroups:
            nrgygrps = nrgygrps+line+' '
    if (type=='md_single_point'):
        dict = {
                'comment':'; Output control options:',
                'nstxout':'0',
                'nstvout':'0',
                'nstfout':'0',
                'nstlog':'0',
                'nstenergy':'0',
                'nstxtcout':'0',
                'xtc_precision':'1000',
                'energygrps':nrgygrps
                }
    if (type=='md_move'):
        dict = {
                'comment':'; Output control options:',
                'nstxout':'0',
                'nstvout':'0',
                'nstfout':'0',
                'nstlog':'0',
                'nstenergy':'0',
                'nstxtcout':'0',
                'xtc_precision':'1000',
                'energygrps':nrgygrps
                }
    if (type=='hybrid' or type=='hybrid_sp'):
        dict = {
                'comment':'; Output control options:',
                'nstxout':str(nsteps),
                'nstvout':str(nsteps),
                'nstfout':str(nsteps),
                'nstlog':'0',
                'nstenergy':str(nsteps),
                'nstxtcout':str(nsteps),
                'xtc_precision':'1000',
                'energygrps':nrgygrps
                }
    if (type=='md'):
        freq = nsteps/10
        dict = {
                'comment':'; Output control options:',
                'nstxout':str(freq),
                'nstvout':str(freq),
                'nstfout':str(freq),
                'nstlog':str(freq),
                'nstenergy':str(freq),
                'nstxtcout':str(freq),
                'xtc_precision':'1000',
                'energygrps':'System'
                }
    return dict

def default_nsearch_dict_mdp():
    dict ={
           'comment':'; Neighborsearching parameters:',
           'rlist':'1.2'
           }
    return dict

def default_coul_vdw_dict_mdp():
    dict = {
            'comment':'; Options for electrostatics and vdw:',
            'coulombtype':'shift',
            'rcoulomb'  :'1.2',
            'epsilon_r' :'15',
            'vdw_type'  :'shift',
            'rvdw_switch':'0.9',
            'rvdw'      :'1.2',
            'fourier_nx':'10',
            'fourier_ny':'10',
            'fourier_nz':'10'
            }
    return dict

def deafult_weak_coupling_dict_mdp(type,ensemble,pressure,temperature):
    dict = {}
    if (type=='em' or type=='em_posre'):
        dict = {
                'comment':'; Options for weak coupling algorithms:',
                'tcoupl':'no',
                'Pcoupl':'no'
                }
    if (type=='md_posre' or type=='md_pull_eq' or type=='md_pull_prod' or type=='md_umbr_prod'):
        dict = {
                'comment':'; Options for weak coupling algorithms:',
                'tcoupl':'Berendsen',
                'tc-grps':'Protein Non-Protein',
                'tau_t':'0.1 0.1',
                'ref_t':'303 303',
                'Pcoupl':'Berendsen',
                'Pcoupltype':'isotropic',
                'tau_p':'0.5',
                'compressibility':'1e-5',
                'ref_p':'300.0',
                'andersen_seed':'815131'
                }
    if (type=='md' or
        type=='md_single_point' or
        type=='md_move' or
        type=='hybrid' or
        type=='hybrid_sp'):
        if(ensemble=='NVT'):
            pcoupl = 'no'
            press  = str(pressure)
            temp   = str(temperature)
        if(ensemble=='NPT'):
            pcoupl = 'Berendsen'
            press  = str(pressure)
            temp   = str(temperature)
        dict = {
                'comment':'; Options for weak coupling algorithms:',
                'tcoupl':'Berendsen',
                'tc-grps':'System',
                'tau_t':'0.1',
                'ref_t':temp,
                'Pcoupl':pcoupl,
                'Pcoupltype':'isotropic',
                'tau_p':'0.5',
                'compressibility':'1e-5',
                'ref_p':press,
                'andersen_seed':'815131'
                }
    return dict

def default_constraints_dict_mdp():
    dict = {}
    if (type=='md_pull_eq' or type=='md_pull_prod' or type=='md_umbr_prod'):
        dict = {
                'comment':'; Constraints:',
                'constraints':'none',
                'constraint_algorithm':'lincs',
                'lincs_iter':'4',
                'unconstrained-start':'no'
                }
    else:
        dict = {
                'comment':'; Constraints:',
                'constraints':'none',
                'constraint_algorithm':'lincs',
                'lincs_iter':'4'
                }
    return dict

def default_genvel_dict_mdp(type,temperature):
    dict = {}
    if (type=='md_posre' or type=='md_pull_eq'):
        dict = {
                'comment':'; Generate velocities for startup run:',
                'gen_vel':'yes',
                'gen_temp':'303',
                'gen_seed':'171533'
                }
    if (type=='md_pull_prod' or type=='md_umbr_prod'):
        dict = {
                'comment':'; Generate velocities for startup run:',
                'gen_vel':'no',
                'gen_temp':'303',
                'gen_seed':'-1'
                }
    temp = str(temperature)
    if (type=='md'):
        dict = {
                'comment':'; Generate velocities for startup run:',
                'gen_vel':'no',
                'gen_temp':temp,
                'gen_seed':'171533'
                }
    if (type=='md_single_point'):
        dict = {
                'comment':'; Generate velocities for startup run:',
                'gen_vel':'no',
                'gen_temp':temp,
                'gen_seed':'171533'
                }
    if (type=='md_move'):
        dict = {
                'comment':'; Generate velocities for startup run:',
                'gen_vel':'no',
                'gen_temp':temp,
                'gen_seed':'171533'
                }
    if (type=='hyrbid' or type=='hyrbid_sp'):
        dict = {
                'comment':'; Generate velocities for startup run:',
                'gen_vel':'no',
                'gen_temp':temp,
                'gen_seed':'171533'
                }
    return dict

def default_pull_dict_mdp(type,dcom='',groups=['A','B']):
    dict = {}
    if (type=='md_pull_eq' or type=='md_pull_prod'):
        dict = {
                'comment':'; Pull stuff:',
                'pull':'constraint',
                'pull_geometry':'distance',
                'pull_group0':groups[0],
                'pull_group1':groups[1],
                'pull_nstxout':'0',
                'pull_nstfout':'5',
                'pull_init1':dcom
                }
    if (type=='md_umbr_prod'):
        dict = {
                'comment':'; Pull stuff:',
                'pull':'umbrella',
                'pull_geometry':'distance',
                'pull_dim':'Y Y Y',
                'pull_start':'yes',  # define initial COM distance > 0
                'pull_ngroups':'1',
                'pull_nstxout':'5',
                'pull_nstfout':'5',
                'pull_group0':groups[0],
                'pull_group1':groups[1],
                'pull_rate1':'0.00', # 0.0 nm per ps
                'pull_k1':'1000'    # kJ mol^-1 nm^-2
                }
    return dict

def write_mdp_dict(f,dict):
    if ('comment' in dict):  # Modified from (dict.has_key('comment'))
        f.write(dict['comment']+'\n')
        del dict['comment']
    if(len(dict)==0):
        return
    onetab   = '\t'
    twotab   = '\t\t'
    threetab = '\t\t\t'
    for key,value in dict.items():  # Changed iteritems() to items()
        tab  = threetab
        if (len(key)>7):
            tab = twotab
        elif (len(key)>15):
            tab = onetab
        f.write(key+tab+'= '+value+'\n')
    f.write('\n')

def write_mdp_header(f,mdpfile):
    f.write(';\n')
    f.write(';'+' '+mdpfile+'\n')
    f.write(';'+' Generated automatically by: '+os.path.abspath(__file__)+'\n')
    user = os.path.abspath(__file__).strip().split('/')[2]
    f.write(';'+' User: '+user+'\n')
    now = dt.datetime.now()
    f.write(';'+' Date: '+str(now.date())+'\n')
    f.write(';'+' Time: '+str(now.time())+'\n')
    f.write(';\n\n')

def generate_mdp(filebase,
                 dcom='',
                 groups=['A','B'],
                 nsteps=0,
                 pressure=1.0,
                 temperature=303,
                 ensemble='NVT',
                 extragroups=[]):
    ext = 'mdp'
    mdpfile = filebase+'.'+ext
    f = open(mdpfile,'w')
    write_mdp_header(f,mdpfile)
    preprocess_dict = default_preprocess_dict_mdp(filebase)
    write_mdp_dict(f,preprocess_dict)
    run_control_dict = default_run_control_dict_mdp(filebase,nsteps)
    try:
        if run_control_dict['integrator']=='md':
            print("Writing mdp for", run_control_dict['nsteps'], 'steps at', run_control_dict['dt'], 'ps timestep; total', float(run_control_dict['nsteps'])*float(run_control_dict['dt'])/1e3, 'ns runtime')
    except KeyError:
        # silently ignore
        pass
    write_mdp_dict(f,run_control_dict)
    output_ctrl_dict = default_output_ctrl_dict_mdp(filebase,nsteps,extragroups)
    write_mdp_dict(f,output_ctrl_dict)
    nsearch_dict = default_nsearch_dict_mdp()
    write_mdp_dict(f,nsearch_dict)
    coul_vdw_dict = default_coul_vdw_dict_mdp()
    write_mdp_dict(f,coul_vdw_dict)
    weak_coupling_dict = deafult_weak_coupling_dict_mdp(filebase,ensemble,pressure,temperature)
    write_mdp_dict(f,weak_coupling_dict)
    genvel_dict = default_genvel_dict_mdp(filebase,temperature)
    write_mdp_dict(f,genvel_dict)
    pull_dict = default_pull_dict_mdp(filebase,dcom,groups)
    write_mdp_dict(f,pull_dict)
    f.close()
    return mdpfile
