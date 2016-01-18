#!/usr/bin/env python
import os,sys,glob,shlex,optparse,subprocess


def DMRGinit():
    usage = """%prog [options] INPUT 
Initialize a DMRG calculation with the CI-DEAS procedure starting from the given input file.
Use --launcher="mpirun OPTIONS" with the appropriate OPTIONS to run an mpi job.
"""
    parser = optparse.OptionParser(usage)

    parser.add_option('-v','--verbose', action='store_true', dest="verbose",
		     help='Print more information',default=False)

    group = optparse.OptionGroup(parser, 'General options')
    group.add_option('--launcher', type='string', default=None,dest="launcher",
		     help='Command used to launch the main executable, i.e. "OMP_NUM_THREADS=2"',metavar="COMMAND")
    group.add_option('--executable', type='string', default=None, dest="executable",
                     help='Path to the main executable',metavar="FILE")
    parser.add_option_group(group)

    group = optparse.OptionGroup(parser, 'File handling  ')
    group.add_option('--g','--get', type='string', default=None,dest="getstring",
		     help='command with no action (only available for compatibility reasons)',metavar="COMMAND")
    group.add_option('--i','--put', type='string', default=None,dest="putstring",
		     help='command with no action (only available for compatibility reasons)',metavar="COMMAND")
    parser.add_option_group(group)

    # Read .rc files
    rcargs = []
    rcpaths = ['dmrgrc','dmrgrc',
               os.path.join(os.path.expanduser(''),'dmrgrc')]
    for path in rcpaths:
        try:
            f = open(path,'r')
            for l in f.readlines():
                l = l.strip()
                if len(l) > 0:
                    if l[0] != '#':
                        for a in shlex.split(l):
                            rcargs.append(os.path.expandvars(a))
            break
        except Exception:
            pass

    (options, args) = parser.parse_args(rcargs+sys.argv[1:])

    return (options.verbose,options.launcher)


def write_di(settings, di):
   difile = open(settings['difile'], 'w')
   fread = open(di,'r')
   for line in fread:
      if line.strip():
         words = line.split()
         if not words[0] in settings:
            difile.write(line)
         else:
            difile.write(words[0] + ' = ' + str(settings[words[0]])+ '\n')
   fread.close()
   difile.close()
   return

def set_fullexe(exename,launcher):
    if launcher is None:
        return (exename+' ')
    else:
        return (' ' + launcher + ' ' + exename + ' ')

#The actual work starts here

(verbose,launcher) = DMRGinit()

nargs = len(sys.argv)

settings = {}
diopen = open(sys.argv[nargs-1],'r')
for line in diopen:
   if line.strip(): #check for empty lines
      words = line.split()
      if words[0] == 'chkpfile':
         settings['chkpfile'] = words[2]
      if words[0] == 'resultfile':
         settings['resultfile'] = words[2]
      if words[0] == 'symmetry':
         settings['symmetry'] = words[2]
diopen.close()

settings['max_bond_dimension'] = 200
settings['nsweeps'] = 2
settings['difile'] = 'di_init'
settings['old_di'] = sys.argv[nargs-1]


if settings['symmetry'][1:-1] == 'su2u1pg':
   settings['det2mps'] = 'det2mps_su2u1pg'
elif settings['symmetry'][1:-1] == 'su2u1':
   settings['det2mps'] = 'det2mps_su2u1'
elif settings['symmetry'][1:-1] == '2u1pg':
   settings['det2mps'] = 'det2mps_2u1pg'
elif settings['symmetry'][1:-1] == '2u1':
   settings['det2mps'] = 'det2mps_2u1'


write_di(settings, sys.argv[nargs-1])

try:
    return_calc = subprocess.Popen([set_fullexe("dmrg ",launcher) +" "+settings['difile']],shell=True)
    return_calc.wait()
    try:
        return_meas = subprocess.Popen([set_fullexe("dmrg_meas ",launcher)+settings['difile']],shell=True)
        return_meas.wait()
        try:
            return_det2mps = subprocess.Popen([settings['det2mps']+' '+settings['old_di']],shell= True)
            return_det2mps.wait()
        except subprocess.CalledProcessError as e:
            print ' CI-DEAS MPS instantiation failed with ', e
    except subprocess.CalledProcessError as e:
        print 'DMRG measurements failed with ', e
except subprocess.CalledProcessError as e:
    print 'DMRG failed with ', e
