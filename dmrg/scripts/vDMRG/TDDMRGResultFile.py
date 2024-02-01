#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#/**
# * @file
# * @copyright This code is licensed under the 3-clause BSD license.
# *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
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

class VibrationalCalculationTypeError(Exception):
    """Raised when the DMRG calculation is not a vDMRG one"""
    pass

class inputError(Exception):
    """Raised if there are incoherence in the data"""
    pass

class TDDMRGMeasurement(object):

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
            if not vibrationalTypeName == 'watson':
                sys.exit("Only Watson-based TD-DMRG calculations supported")
            # Extracts the key parameters
            self.overallSize = self.h5pyfile['parameters']['L'][()]
            self.nMaxParameter = self.h5pyfile['parameters']['Nmax'][()]
            self.nSweeps = self.h5pyfile['parameters']['nsweeps'][()]
        except:
            raise h5Error

    def getNumberOfModes(self):
        """Getter for the number of modal bases"""
        return self.overallSize

    def getLatticeSize(self):
        """Getter for the overall number of sites"""
        return self.overallSize

    def getNumberOfSweeps(self):
        """Getter for the number of sweeps"""
        return self.nSweeps

    def extractModeExcitationDegree(self, nElements: int):
        """Extracts the one-mode RDM matrix"""
        excitationDegree = []
        numberOfSweeps = nElements
        if numberOfSweeps == 0:
            numberOfSweeps = self.nSweeps
        # Loop over the number of modes
        for iMode in range(self.overallSize):
            populationDynamics = []
            for iStep in range(numberOfSweeps):
                populationDynamics.append(self.h5pyfile['spectrum']['iteration'][str(iStep)]['results']['ExcitationMode'+str(iMode)]['mean']['value'][0][0][0])
            excitationDegree.append(populationDynamics)
        return excitationDegree

