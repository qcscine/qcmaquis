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

