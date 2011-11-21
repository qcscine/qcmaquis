#!/usr/bin/env python

import sys, os

p = os.popen('svnversion -n '+sys.argv[1])
rev = p.readline()
p.close()

print sys.argv

fin = open(sys.argv[2], 'r')
fout = open(sys.argv[3], 'w')
for line in fin.readlines():
	if '#cmakedefine DMRG_VERSION "@DMRG_VERSION@"' in line:
		line = '#define DMRG_VERSION "0.1-'+rev+'"\n'
	if '#cmakedefine DMRG_VERSION_MAJOR @DMRG_VERSION_MAJOR@' in line:
		line = '#define DMRG_VERSION_MAJOR 0\n'
	if '#cmakedefine DMRG_VERSION_MINOR @DMRG_VERSION_MINOR@' in line:
		line = '#define DMRG_VERSION_MINOR 1\n'
	if '#cmakedefine DMRG_VERSION_BUILD @DMRG_VERSION_BUILD@' in line:
		line = '#define DMRG_VERSION_BUILD "'+rev+'"\n'
	if '#cmakedefine DMRG_VERSION_STRING "@DMRG_VERSION_STRING@"' in line:
		line = '#define DMRG_VERSION_STRING "DMRG version 0.1-'+rev+'"\n'
	fout.write(line)
fin.close()
fout.close()
