#!/usr/bin/env python

#######################################################
#
#Driver to automate excited state calculations
#from templated input files and varying number
#of excited states
#
#Christopher J. Stein, ETH Zurich, October 2016
######################################################


import os
import jinja2
import sys
import subprocess
import numpy as np
import pydmrg


#Declaration of variables

num_ex_states = int(sys.argv[1])+1

integral_file = '"/home/eth/steinc/small_projects/watson/excited/CH2S_VibHam"'
maquis_path = '/home/eth/steinc/build/maquis_watson/applications/dmrg/dmrg'
input_file = 'di'


L                  = 6
nsweeps            = 5
max_bond_dimension = 10
Nmax               = 9
ortho_states       = '""'



for state in range(1,num_ex_states,1):
   resultfile = '"result_'+str(state)+'.h5"'
   chkpfile = '"chkp_'+str(state)+'.h5"'
   n_ortho_states = state-1

   #write QCMaquis input file from template
   templateLoader = jinja2.FileSystemLoader(searchpath=os.path.dirname(os.path.realpath(__file__)))
   templateEnv = jinja2.Environment(loader=templateLoader)
   TEMPLATE_FILE = 'Templates/dmrg.input'
   t = templateEnv.get_template( TEMPLATE_FILE )

   dmrg_input = open(input_file,'w')
   dmrg_input.write(t.render(L = L, Nmax = Nmax, nsweeps = nsweeps, max_bond_dimension = max_bond_dimension,
                             integral_file = integral_file, resultfile = resultfile, chkpfile = chkpfile,
                             n_ortho_states = n_ortho_states, ortho_states = ortho_states))
   dmrg_input.close()

   #update string of orthogonal states
   if n_ortho_states == 0:
      ortho_states = ortho_states[:-1]+chkpfile[1:]
   else:
      ortho_states = ortho_states[:-1]+','+chkpfile[1:]

   #run QCMaquis DMRG calculation
   try:
      calc = subprocess.Popen([maquis_path + " " + input_file],shell=True)
      calc.wait()
   except subprocess.CalledProcessError as e:
      print 'DMRG calculation failed with ', e

#print energies to screen and file
outfile = open('out','w')

for state in range(1,num_ex_states,1):
   resultfile = 'result_'+str(state)+'.h5'
   energy_data = pydmrg.LoadDMRGSweeps([resultfile],['Energy'])
   sweeps = []
   for sw in energy_data[0]:
      sweeps += list(sw[0].y)
   ydata = np.array(sweeps)

   min_energy = np.min(ydata.real)
   min_energy_last_sweep = np.min(ydata[(2*(L-1))*(nsweeps-1)-1:(2*(L-1))*nsweeps-1])
   min_energy_second_to_last_sweep = np.min(ydata[(2*(L-1))*(nsweeps-2)-1:(2*(L-1))*(nsweeps-1)-1])

   print "State: ", state
   print "Minimum energy:", min_energy
   print "Minimum energy of last sweep:", min_energy_last_sweep
   if (min_energy_second_to_last_sweep-min_energy_last_sweep) > 0.01:
      print "     CAUTION: energy decreased by more than 0.01 cm^-1 in the last sweep"
      print "     ", min_energy_second_to_last_sweep, "cm^-1  vs. ", min_energy_last_sweep, "cm^-1 "

   outfile.write("State: "+str(state)+"\n")
   outfile.write("Minimum energy: " +str(min_energy)+"\n")
   outfile.write("Minimum energy of last sweep: "+str(min_energy_last_sweep)+"\n")
   if (min_energy_second_to_last_sweep-min_energy_last_sweep) > 0.01:
      outfile.write("     CAUTION: energy decreased by more than 0.01 cm^-1 in the last sweep!\n")
      outfile.write("     "+str(min_energy_second_to_last_sweep)+"cm^-1  vs. "+str(min_energy_last_sweep)+"cm^-1\n\n\n")
   else:
      outfile.write("\n\n")

outfile.close()


