#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#/**
# * @file
# * @copyright This code is licensed under the 3-clause BSD license.
# *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
# *            See LICENSE.txt for details.
# */

from enum import Enum
import h5py
import math
import numpy as np
import sys
from ExcitonicTDResultsFile import ResultFileVibrationaleMeasurement

class ResultsFilePopulationMeasurement(ResultFileVibrationaleMeasurement):
    
    def extractNumMolecules(self):
        return self.h5pyfile()['parameters']['vibronic_num_molecules'][()]

    def extractNumSweeps(self):
        return self.h5pyfile['parameters']['nsweeps'][()]-1
    
    def extractNumModes(self):
         return self.h5pyfile['parameters']['vibronic_num_vibmodes'][()]

    def extractPopulations(self):
        nsweeps = self.extractNumSweeps()
        n_elestates = self.extractNumMolecules-1
        pop = np.zeros((self.num_elestates+1, nsweeps), dtype=float) #stores populations for all electornic states
        pop_error = np.zeros((self.num_elestates+1, nsweeps), dtype=float) #stores population error for all electornic states
        for iteration_idx in range(nsweeps): #group index
            for ele_idx in range(self.num_elestates+1): #loop over electronic states   
                a = self.h5pyfile['spectrum']['iteration'][str(iteration_idx)]['results']['PopulationState'+str(ele_idx)]['mean']['value'][0]
                pop[ele_idx, iteration_idx] = a[0,0]
                pop_error[ele_idx, iteration_idx] = a[0,1]
        return pop, pop_error

    def extractAutocorrelation(self):
        nsweeps = self.extractNumSweeps()
        autocorrelation = np.zeros((nsweeps,2), dtype=float)

        for iteration_idx in range(nsweeps): #group index
                entry = self.h5pyfile['spectrum']['iteration'][str(iteration_idx)]['results']['Autocorrelation']['mean']['value'][0]
                autocorrelation[iteration_idx] = entry
        return autocorrelation #first entry corresponds to the real part and the second entry to the imaginary part
    

    def extractDisplacements(self):
         nsweeps = self.extractNumSweeps()
         num_monomers = self.h5pyfile['parameters']['vibronic_num_molecules'][()]
         num_modes = self.h5pyfile['parameters']['vibronic_num_vibmodes'][()]
         num_connectingmodes = self.h5pyfile['parameters']['vibronic_num_connectingmodes'][()]
         displacement_local = np.zeros((nsweeps, num_monomers, num_modes-num_connectingmodes,2), dtype=float)
         displacement_connecting = np.zeros((nsweeps, num_monomers-1, num_connectingmodes, 2), dtype=float)

         for iteration_idx in range(nsweeps): #group index
                for monomer_idx in range(num_monomers):
                     for mode_idx in range(num_modes):
                        if ((monomer_idx == num_monomers-1) and (mode_idx >= (num_modes-num_connectingmodes))):
                            break
                        entry = self.h5pyfile['spectrum']['iteration'][str(iteration_idx)]['results']['Displacement'+str(monomer_idx)+'Mode'+str(mode_idx)]['mean']['value'][0]
                        if (mode_idx < (num_modes-num_connectingmodes)):
                            displacement_local[iteration_idx][monomer_idx][mode_idx] = entry
                        else:
                            displacement_connecting[iteration_idx][monomer_idx][monomer_idx-(num_modes-num_connectingmodes)] = entry
         return displacement_local, displacement_connecting
    
    def extractDisplacement(self):
         nsweeps = self.extractNumSweeps()
         num_monomers = self.h5pyfile['parameters']['vibronic_num_molecules'][()]
         num_modes = self.h5pyfile['parameters']['vibronic_num_vibmodes'][()]
         num_connectingmodes = self.h5pyfile['parameters']['vibronic_num_connectingmodes'][()]
         displacement = np.zeros((nsweeps, num_modes*num_monomers-num_connectingmodes,2), dtype=float)
         count = 0
         switch = 0
         for iteration_idx in range(nsweeps): #group index
                switch = 0
                for monomer_idx in range(num_monomers):
                     for mode_idx in range(num_modes):
                        if ((monomer_idx == num_monomers-1) and (mode_idx >= (num_modes-num_connectingmodes))):
                            switch = 1
                            break
                        entry = self.h5pyfile['spectrum']['iteration'][str(iteration_idx)]['results']['Displacement'+str(monomer_idx)+'Mode'+str(mode_idx)]['mean']['value'][0]
                        displacement[iteration_idx][count] = entry
                        count += 1
                     if(switch==1):
                          break
                count = 0
         return displacement
                        
    def extractDisplacementSquared(self):
         nsweeps = self.extractNumSweeps()
         num_monomers = self.h5pyfile['parameters']['vibronic_num_molecules'][()]
         num_modes = self.h5pyfile['parameters']['vibronic_num_vibmodes'][()]
         num_connectingmodes = self.h5pyfile['parameters']['vibronic_num_connectingmodes'][()]
         displacementSquared = np.zeros((nsweeps, num_modes*num_monomers-num_connectingmodes,2), dtype=float)
         count = 0
         switch = 0
         for iteration_idx in range(nsweeps): #group index
                switch = 0
                for monomer_idx in range(num_monomers):
                     for mode_idx in range(num_modes):
                        if ((monomer_idx == num_monomers-1) and (mode_idx >= (num_modes-num_connectingmodes))):
                            switch = 1
                            break
                        entry = self.h5pyfile['spectrum']['iteration'][str(iteration_idx)]['results']['DisplacementSquared'+str(monomer_idx)+'Mode'+str(mode_idx)]['mean']['value'][0]
                        displacementSquared[iteration_idx][count] = entry
                        count += 1
                     if(switch==1):
                          break
                count = 0
         return displacementSquared
    
    def extractBondDimension(self):
        nsweeps = self.extractNumSweeps()
        L = self.getLatticeSize()
        BondDimension = np.zeros((nsweeps, (L-1)*2-2), dtype=float) #for excitonicextended Lattice
        for iteration_idx in range(nsweeps): #group index   
            BondDimension[iteration_idx][:] = self.h5pyfile['spectrum']['iteration'][str(iteration_idx)]['results']['BondDimension']['mean']['value']
        return BondDimension
                


