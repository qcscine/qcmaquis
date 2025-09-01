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

class h5Error(Exception):
    """Raised when there is a problem opening the h5 file"""
    pass

class VibronicCalculationTypeError(Exception): #VAL:modified
    """Raised when the DMRG calculation is not a vDMRG one"""
    pass

class inputError(Exception):
    """Raised if there are incoherence in the data"""
    pass

class VibronicCalculationType(Enum): #VAL:modified
    """Enum class representing the type of vibrational calculation"""
    EXCITONIC = 1
    VIBRONIC = 2
    EXCITONICEXTENDED = 3  

class ResultFileVibrationaleMeasurement(object):

    ''' Static member '''
    atol = 1.0E-12

    '''
    Class constructor
    '''
    def __init__(self, name): #VAL:modified
        self.nameFile = name
        # Opens the h5 file. If not possible, raises an error
        
        try:
            self.h5pyfile = h5py.File(self.nameFile)
            # Checks which type of vibrational calculation
            vibronicTypeName = self.h5pyfile['parameters']['MODEL'][()].decode('ascii')
            if (vibronicTypeName == 'excitonic'):
                self.__vibronicCalculationType = VibronicCalculationType.EXCITONIC
            elif (vibronicTypeName == 'vibronic'):
                self.__vibronicCalculationType = VibronicCalculationType.VIBRONIC
            elif (vibronicTypeName == 'excitonicextended'):
                self.__vibronicCalculationType = VibronicCalculationType.EXCITONICEXTENDED
            else:
                raise VibronicCalculationTypeError
            # Extracts the key parameters
            self.overallSize = self.h5pyfile['parameters']['L'][()]
            nMaxParameter = self.h5pyfile['parameters']['Nmax'][()]
            self.numBasisFunctions = [nMaxParameter]*self.overallSize
        except:
            raise h5Error
            
    def getNumberOfModes(self):
        """Getter for the number of modal bases"""
        return len(self.numBasisFunctions)

    def getNumberOfModals(self):
        """Getter for the number of modal bases"""
        return self.numBasisFunctions

    def getLatticeSize(self):
        """Getter for the overall number of sites"""
        return self.overallSize

    def __extractOneModeRDM(self):
        """Extracts the one-mode RDM matrix"""
        oneModeRDM = np.zeros((self.overallSize, self.overallSize))
        ihd5 = self.h5pyfile['spectrum']['results']['onemodeRDM']
        for idx, iLabel in enumerate(ihd5['labels_num']):
            oneModeRDM[iLabel[0], iLabel[1]] = abs(ihd5['mean']['value'][0, idx])
        print(oneModeRDM[0,0])
        print(oneModeRDM[1,0])
        print(oneModeRDM[0,1])
        print(oneModeRDM[1,1])
        return oneModeRDM

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

    def __checkAndExtractTwoParticleQuantity(resultFile, label):
        for idx, i in enumerate(resultFile['labels_num']):
            if i[0] == label[0] and i[1] == label[1]:
                return resultFile['mean']['value'][0,idx]
        return None

    def __extractTwoModalEntropy(self):
        """Extracts the two-modals entropy for each modal"""
        twoModalEntropy = np.zeros((self.overallSize, self.overallSize))
        twoModalEntropy00 = np.zeros((self.overallSize, self.overallSize))
        twoModalEntropy11 = np.zeros((self.overallSize, self.overallSize))
        twoModalEntropy22 = np.zeros((self.overallSize, self.overallSize))
        twoModalEntropy33 = np.zeros((self.overallSize, self.overallSize))
        twoModalEntropy12 = np.zeros((self.overallSize, self.overallSize))
        twoModalEntropy21 = np.zeros((self.overallSize, self.overallSize))
        #
        for iMat, ihd5 in zip([twoModalEntropy00, twoModalEntropy11, twoModalEntropy22, twoModalEntropy33, twoModalEntropy12, twoModalEntropy21],
                              [self.h5pyfile['spectrum']['results']['twomodeRDM_00'], self.h5pyfile['spectrum']['results']['twomodeRDM_11'],
                               self.h5pyfile['spectrum']['results']['twomodeRDM_22'], self.h5pyfile['spectrum']['results']['twomodeRDM_33'],
                               self.h5pyfile['spectrum']['results']['twomodeRDM_12'], self.h5pyfile['spectrum']['results']['twomodeRDM_21']]):
            for idx, iLabel in enumerate(ihd5['labels_num']):
                iMat[iLabel[0], iLabel[1]] = ihd5['mean']['value'][0, idx]

        for i in range(self.overallSize):
            for j in range(self.overallSize):
                eigenvalue0 = twoModalEntropy00[i, j]
                eigenvalue1 = (twoModalEntropy11[i, j] + twoModalEntropy22[i, j] + np.sqrt((twoModalEntropy11[i, j] - twoModalEntropy22[i, j])**2 + 4.*twoModalEntropy21[i, j]*twoModalEntropy12[i, j]))/2.
                eigenvalue2 = (twoModalEntropy11[i, j] + twoModalEntropy22[i, j] - np.sqrt((twoModalEntropy11[i, j] - twoModalEntropy22[i, j])**2 + 4.*twoModalEntropy21[i, j]*twoModalEntropy12[i, j]))/2.
                eigenvalue3 = twoModalEntropy33[i, j]
                for iEigen in [eigenvalue0, eigenvalue1, eigenvalue2, eigenvalue3]:
                    if abs(iEigen) > ResultFileVibrationaleMeasurement.atol:
                        twoModalEntropy[i, j] += -iEigen*math.log(iEigen)
        return twoModalEntropy

    def getOneModeRDM(self):
        """
        Gets the one-modal entropy as a "full" list
        """
        return self.__extractOneModeRDM()

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
        """
        Calculates the mutual information
        """
        mutual_information = np.zeros((self.overallSize, self.overallSize))
        singleEntropy = self.__extractOneModalEntropy()
        twoEntropy = self.__extractTwoModalEntropy()
        for i in range(self.overallSize):
            for j in range(i+1, self.overallSize):
                mutual_information[i][j] = 0.5*(singleEntropy[i] + singleEntropy[j] - twoEntropy[i, j])
                mutual_information[j][i] = 0.5*(singleEntropy[i] + singleEntropy[j] - twoEntropy[i, j])
        return mutual_information
