#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2022 Laboratory for Physical Chemistry, ETH Zurich
#*               2022- by Alberto Baiardi <abaiardi@ethz.ch>
#*
#* This software is part of the ALPS Applications, published under the ALPS
#* Application License; you can use, redistribute it and/or modify it under
#* the terms of the license, either version 1 or (at your option) any later
#* version.
#*
#* You should have received a copy of the ALPS Application License along with
#* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
#* available from http://alps.comp-phys.org/.
#*
#* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
#* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
#* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
#* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#* DEALINGS IN THE SOFTWARE.
#*
#*****************************************************************************

from enum import Enum
import h5py
import math
import numpy as np
import sys

class h5Error(Exception):
    """Raised when there is a problem opening the h5 file"""
    pass

class VibrationalCalculationTypeError(Exception):
    """Raised when the DMRG calculation is not a vDMRG one"""
    pass

class inputError(Exception):
    """Raised if there are incoherence in the data"""
    pass

class VibrationalCalculationType(Enum):
    """Enum class representing the type of vibrational calculation"""
    CANONICAL = 1
    NMODE = 2

class ResultFileVibrationaleMeasurement(object):

    ''' Static member '''
    atol = 1.0E-12

    '''
    Class constructor
    '''
    def __init__(self, name):
        self.nameFile = name
        # Opens the h5 file. If not possible, raises an error
        try:
            self.h5pyfile = h5py.File(self.nameFile)
            # Checks which type of vibrational calculation
            vibrationalTypeName = self.h5pyfile['parameters']['MODEL'][()].decode('ascii')
            if (vibrationalTypeName == 'nmode'):
                self.__vibrationalCalculationType = VibrationalCalculationType.NMODE
            elif (vibrationalTypeName == 'watson'):
                self.__vibrationalCalculationType = VibrationalCalculationType.CANONICAL
            else:
                raise VibrationalCalculationTypeError
            # Extracts the key parameters
            self.overallSize = self.h5pyfile['parameters']['L'][()]
            # For the n-mode case, the code needs to extract also 
            if (self.__vibrationalCalculationType == VibrationalCalculationType.NMODE):
                tmpStr = str(self.h5pyfile['parameters']['nmode_num_basis'][()] )[2:-1]
                tmpLst = tmpStr.split(',')
                self.numBasisFunctions = []
                overallSizeCheck = 0
                for i in tmpLst:
                    self.numBasisFunctions.append(int(i))
                    overallSizeCheck += int(i)
                if (overallSizeCheck != self.overallSize):
                    raise inputError
            elif (self.__vibrationalCalculationType == VibrationalCalculationType.CANONICAL):
                nMaxParameter = self.h5pyfile['parameters']['Nmax'][()]
                self.numBasisFunctions = [nMaxParameter]*self.overallSize
        except:
            raise h5Error
    
    def getNumberOfModals(self):
        """Getter for the number of modal bases"""
        return self.numBasisFunctions

    def getLatticeSize(self):
        """Getter for the overall number of sites"""
        return self.overallSize

    def __extractOneModalEntropy(self):
        """Extracts the one-modal entropy for each modal"""
        oneModalEntropy = np.zeros(self.overallSize)
        oneModalMatrix00 = np.zeros(self.overallSize)
        oneModalMatrix11 = np.zeros(self.overallSize)
        # First element of the entanglement matrix
        for iLabel, iElement in zip(self.h5pyfile['spectrum']['results']['onemodalRDM_00']['labels_num'],
                                    self.h5pyfile['spectrum']['results']['onemodalRDM_00']['mean']['value'][0]):
            oneModalMatrix00[iLabel] = iElement
        # Second element of the entanglement matrix
        for iLabel, iElement in zip(self.h5pyfile['spectrum']['results']['onemodalRDM_11']['labels_num'],
                                    self.h5pyfile['spectrum']['results']['onemodalRDM_11']['mean']['value'][0]):
            oneModalMatrix11[iLabel] = iElement
        # Calculates the quantum entropy
        for iCont, (i00, i11) in enumerate(zip(oneModalMatrix00, oneModalMatrix11)):
            if abs(i00) > ResultFileVibrationaleMeasurement.atol:
                oneModalEntropy[iCont] += -i00*math.log(i00)
            if abs(i11) > ResultFileVibrationaleMeasurement.atol:
                oneModalEntropy[iCont] += -i11*math.log(i11)
        return oneModalEntropy

    # # Extracts the two-modals entropy for each modal
    # def extractTwoModalEntropy(self):
    #     ret = np.zeros(self.overallSize, self.overallSize)
    #     for jj in range(self.overallSize):
    #         for kk in range(jj+1, self.overallSize):
    #             eigenvalues = np.zeros((4))
    #             eigenvalues[0] = self.h5pyfile['spectrum']['results']['twomodeRDM_'+str(jj)+"_"+str(kk)+"_00"]['mean']['value'][0]
    #             eigenvalues[3] = self.h5pyfile['spectrum']['results']['twomodeRDM_'+str(jj)+"_"+str(kk)+"_33"]['mean']['value'][0]
    #             a = self.h5pyfile['spectrum']['results']['twomodeRDM_'+str(jj)+"_"+str(kk)+"_11"]['mean']['value'][0]
    #             b = self.h5pyfile['spectrum']['results']['twomodeRDM_'+str(jj)+"_"+str(kk)+"_12"]['mean']['value'][0]
    #             c = self.h5pyfile['spectrum']['results']['twomodeRDM_'+str(jj)+"_"+str(kk)+"_21"]['mean']['value'][0]
    #             d = self.h5pyfile['spectrum']['results']['twomodeRDM_'+str(jj)+"_"+str(kk)+"_22"]['mean']['value'][0]
    #             eigenvalues[1] = (a + d + np.sqrt((a - d)**2 + 4.0*b*c))/2.0
    #             eigenvalues[2] = (a + d - np.sqrt((a - d)**2 + 4.0*b*c))/2.2
    #             ret[jj][kk] = self.__calculateEntropy(eigenvalues)
    #             ret[kk][jj] = ret[jj][kk]
    #     return ret
    
    def getOneModalEntropy(self):
        """
        Gets the one-modal entropy as a "full" list
        """
        return self.__extractOneModalEntropy()

    def getSeparatedOneModalEntropy(self):
        """
        Gets the one-modal entropy as a list of lists, one per mode
        """
        lstOfList = []
        for i in range(len(self.numBasisFunctions)):
            lstOfList.append([])
        lstOfEntropies = self.__extractOneModalEntropy()
        iMode = 0
        for i in range(self.overallSize):
            lstOfList[iMode].append(lstOfEntropies[i])
            if (len(lstOfList[iMode]) >= self.numBasisFunctions[iMode]):
                iMode = iMode + 1
        return lstOfList

    def getMutualInformation(self):
        """Calculates the mutual information"""
        mutual_information = np.zeros((self.overallSize, self.overallSize))
        single_entropy = self.__extractOneModalEntropy()
        # two_entropy = self.extractTwoModalEntropy()
        # for i in range(self.overallSize):
        #     for j in range(i+1, self.overallSize):
        #         result = self.__calculateMutualInformation(single_entropy, two_entropy, i, j)
        #         mutual_information[i][j] = result
        #         mutual_information[j][i] = result
        return mutual_information

    # def __calculateEntropy(self, w):
    #     w = w[w > ResultFileVibrationaleMeasurement.atol]
    #     lnw = np.log(w)
    #     return -1.0*np.sum(w*lnw)

    # def __calculateMutualInformation(single_ent, two_ent, i, j):
    #     if i == j:
    #         return 0
    #     else:
    #         return 0.5*(single_ent[i] + single_ent[j] - two_ent[i][j])
