#!/usr/bin/env python3

#/**
# * @file
# * @copyright This code is licensed under the 3-clause BSD license.
# *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
# *            See LICENSE.txt for details.
# */

""" This module emulates pyalps functions that are required in QCMaquis, but
    is compatible with python 3.
    Since we no longer do not want to depend on pyalps, we needed a replacement.
"""
import h5py, re

class DataSet:
    """Emulates the ALPS DataSet class but is simpler
    """
    def __init__(self,x=None,y=None):
        if x is None: self.x = []
        else:         self.x = x
        if y is None: self.y = numpy.array([])
        else:         self.y = y

def loadEigenstateMeasurements(inputfile, what):
    """ loads a measurement in what from QCMaquis [ALPS] HDF5 file.
        Unlike the real pyalps, we can only load one file at a time here
    """

    # emulate the behaviour of pyalps that expects a list of strings instead of one string
    # so if we get a list, we care only about the first element

    if isinstance(inputfile,list):
        inputfile = inputfile[0]

    if isinstance(what,list):
        what = what[0]

    f=h5py.File(inputfile,'r+')
    data = f.get('/spectrum/results/' + what)

    # Handle the case if the measurement does not exist
    # -> data is None

    if data is None:
        raise RuntimeError("Could not find measurement " + what + " in file " + inputfile + ". Did you perform the measurement?")

    labels = []
    for l in data["labels"]:
        labels.append(list(map(lambda x: int(x),re.findall('[0-9]+',str(l),re.DOTALL))))

    # For compatibility with pyalps, return the value as list of lists with 1 element (woo!)
    return [[DataSet(labels,data["mean/value"][0])]]

def getParameters(inputfile):
   """ Extracts the parameters from HDF5 files
   """
   if isinstance(inputfile, list):
       inputfile = inputfile[0];

   f=h5py.File(inputfile,'r+')

   # read "/parameters" from HDF5 file
   data=f.get("parameters")

   # into a dictionary
   params_dict = {}
   for k in data:
       d = data[k][()]
       if isinstance(d,bytes):
           d=d.decode()
       params_dict[k] = d

   # return a list with one element for compatibility with real pyalps
   return [params_dict]
