#!/usr/bin/env python

# Script that changes the names of QCMaquis checkpoint files inside OpenMOLCAS *.h5 files
# Required for MPSSI calculations

# Copyright 2016-2020 Leon Freitag

import h5py
import sys
import re

# The QCMaquis checkpoint record is stored as a fixed-length string of 256 characters.
HDF_STR_LEN=256

def sort_args():
  for i in range(2,len(sys.argv)-1,2):
    a=(int(sys.argv[i]), sys.argv[i+1])
    yield a

def main():
  quickmode = False
  usage = """
    This script patches the MOLCAS rasscf.h5 file with new names for QCMaquis checkpoint states.
    Use it e.g. for a RASSI calculation in a new project after several different DMRG-SCF
    calculations in separate project folders.

    Usage: qcm_checkpoint_rename.py <rasscf.h5 filename> <state#1> <new name for state #1> <state#2> <new name for state #2>...

    (Warning -- state numbering begins from zero!)

    or

    qcm_checkpoint_rename.py <rasscf.h5 filename> -q

    which enforces "quick" rename -- all checkpoint names in the file will be adjusted to the project
    name of the rasscf.h5 name.
  """

  if len(sys.argv) < 3:
    print(usage)
    sys.exit(0)

  if sys.argv[2] == '-q':
    quickmode = True
# arguments must come in tuples + there's argv[0] -- total length of sys.arg mod 2 must be 1
  elif (len(sys.argv) % 2):
    print("Error: Number of arguments must be even -- pairs of state numbers and names.")
    print(usage)
    sys.exit(1)

  try:
    f = h5py.File(sys.argv[1],'r+')
  except:
    print('Error opening HDF5 file %s' % sys.argv[1])
    sys.exit(1)

  try:
    data = f.get('QCMAQUIS_CHECKPOINT')
  except:
    print('QCMAQUIS_CHECKPOINT record not found. Not a rasscf.h5 file or generated with an old MOLCAS version.')
    sys.exit(1)

  if quickmode:
    name2=sys.argv[1]
    name=sys.argv[1].replace(".rasscf.h5","")
    if (name2 == name):
      name=sys.argv[1].replace(".h5","")
    for i,d in enumerate(data):
      substr = re.sub(b"^.+\.checkpoint_state", name.encode() + b".checkpoint_state", d)
      print('Renaming checkpoint from %s to %s' % (d.decode('utf-8'),substr.decode('utf-8')))
      data[i] = substr.ljust(HDF_STR_LEN)

    for state in sort_args():
      print('Renaming state %i checkpoint %s to %s' % (state[0], data[state[0]], state[1]))
      data[state[0]] = state[1].ljust(HDF_STR_LEN)

  f.close()
if __name__ == "__main__":
  main()
