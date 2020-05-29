"""
OsTools.py

Modified for python 3:
10-12-2019 By Arriën S. Rauh
- Changed print statements for python 3x
- Added filename
"""
import os

Verbose   = True # sets logging of system command


def system(command, options=None, stdin=None, log=None, err=None):
    command = command.replace(' ','\ ')
    if Verbose:
        print(command, options, stdin, log)
    if options:
        command += ' '+' '.join(options)
    if log: # redirect stdout to log file
        command += ' > '+log
    if err: # redirect stderr to err file
        command += ' 2> '+err
    elif log: # if we're redirecting stdin, but not stderr
              # send stderr also to the same log file
        command += ' 2>&1'
    if stdin: # redirect stdin (provide scripted user input)
        command += ' ' + '\n'.join(['<<EOF']+stdin+['EOF'])
    status = os.system(command)
    if Verbose and status:
        print("Failed:", status)
        print("$", command)
    return status


def MakeDir(Dir):
    if (not (os.path.exists(Dir) and os.path.isdir(Dir))):
        os.mkdir(Dir)


def GotoDir(Dir):
    if (not (os.path.exists(Dir) and os.path.isdir(Dir))):
        os.mkdir(Dir)
    os.chdir(Dir)


def cp(source, dest):
    return system('cp \"'+source+'\" \"'+dest+'\"')


def symlink(src,dest):
    if Verbose:
        print("Link", src, dest)
    if os.path.exists(src):
        if os.path.isfile(dest) or os.path.islink(dest):
            os.unlink(dest)
        return os.symlink(src,dest)
    return -1


def checkfile(path,file):
    fullpath = path+'/'+file
    e = os.path.exists(fullpath)
    f = os.path.isfile(fullpath)
    l = os.path.islink(fullpath)
    boOk = e and (f or (f and l))
    return boOk


def put_ifile_in_memory(ifile):
    fi     = open(ifile,"r")
    fi_mem = fi.readlines()
    fi.close()
    length = len(fi_mem)
    for i in range(length):
        l         = fi_mem[i].strip().split()
        fi_mem[i] = l
    return fi_mem


def filename(path_file):
    """
    Extracts filename without extension from full path.
    """
    file = os.path.split(path_file)[-1]
    filename = file.split('.')[0]
    return filename


def file_extension(path_file):
    """
    Extracts file extension from full path.
    """
    file = os.path.split(path_file)[-1]
    extension = file.split('.')[-1]
    return extension
    
    
    
    
    