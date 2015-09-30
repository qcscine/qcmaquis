import sys
import subprocess

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

settings = {}
diopen = open(sys.argv[1],'r')
for line in diopen:
   if line.strip(): #check for empty lines
      words = line.split()
      if words[0] == 'init_chkpfile':
         settings['chkpfile'] = words[2]
      if words[0] == 'init_resultfile':
         settings['resultfile'] = words[2]
      if words[0] == 'symmetry':
         settings['symmetry'] = words[2]
diopen.close()

settings['max_bond_dimension'] = 200
settings['nsweeps'] = 2
settings['difile'] = 'di_init'
settings['old_di'] = sys.argv[1]
settings['res_adapted'] = settings['resultfile'][1:-1]+'_adapted'
settings['chkp_adapted'] = settings['chkpfile'][1:-1]+'_adapted'


if settings['symmetry'][1:-1] == 'su2u1pg':
   settings['transform'] = 'mps_transform_pg'
   settings['new_sym'] = '\'2u1pg\''
   settings['det2mps'] = 'det2mps_su2u1pg'
elif settings['symmetry'][1:-1] == 'su2u1':
   settings['transform'] = 'mps_transform'
   settings['new_sym'] = '\'2u1\''
   settings['det2mps'] = 'det2mps_su2u1'
elif settings['symmetry'][1:-1] == '2u1pg':
   settings['det2mps'] = 'det2mps_2u1pg'
elif settings['symmetry'][1:-1] == '2u1':
   settings['det2mps'] = 'det2mps_2u1'




write_di(settings, sys.argv[1])

if settings['symmetry'][1:-1] == '2u1pg' or settings['symmetry'][1:-1] == '2u1':
   try:
      return_calc = subprocess.Popen(["dmrg "+settings['difile']],shell=True)
      return_calc.wait()
      try:
         return_meas = subprocess.Popen(["dmrg_meas "+settings['difile']],shell=True)
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
else:
   try:
      return_calc = subprocess.Popen(["dmrg "+settings['difile']],shell=True)
      return_calc.wait()
      cp1 = subprocess.Popen(["mv "+settings['chkpfile']+' '+settings['chkp_adapted']], shell=True)
      try:
         print settings['transform']+' '+settings['chkp_adapted']
         return_transform = subprocess.Popen([settings['transform']+' '+settings['chkp_adapted']], shell=True)
         return_transform.wait()
         try:
            settings['symmetry'] = settings['new_sym']
            write_di(settings, sys.argv[1])
            cp2 = subprocess.Popen(["mv "+settings['chkpfile'][1:-4]+'_adapted.*.h5 '+settings['chkpfile'][1:-1]], shell=True)
            return_meas = subprocess.Popen(["dmrg_meas "+settings['difile']],shell=True)
            return_meas.wait()
            try:
               return_det2mps = subprocess.Popen([settings['det2mps']+' '+settings['old_di']],shell= True)
               return_det2mps.wait()
            except subprocess.CalledProcessError as e:
               print ' CI-DEAS MPS instantiation failed with ', e
         except subprocess.CalledProcessError as e:
            print 'DMRG measurements failed with ', e
      except subprocess.CalledProcessError as e:
         print "Transformation from SU2 to 2U1 failed with ", e
   except subprocess.CalledProcessError as e:
      print 'DMRG failed with ', e
