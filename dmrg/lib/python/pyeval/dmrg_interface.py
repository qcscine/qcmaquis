#!/usr/bin/python
import argparse
import math
import numpy as np
import re
import sys
#
# LOCAL VARIABLES
# ---------------
DESCR = '''
Interface from GAUSSIAN anharmonic .log file to Vib-DMRG program
'''
pattern_harmonic       = r'\sMode\(n\)\s+Status\s+E\(harm\)\s+E\(anharm\)' 
pattern_cubic          = r':\s+CUBIC\sFORCE\sCONSTANTS\sIN\sNORMAL\sMODES'
pattern_quartic        = r':\s+QUARTIC\sFORCE\sCONSTANTS\sIN\sNORMAL\sMODES'
pattern_Coriolis       = r':\s+CORIOLIS COUPLINGS'
pattern_harmonic_const = r'''
\s+(?P<index>\d+)\(1\)\s+(?:active|passive)\s+(?P<harm>\d+\.\d+)
'''
pattern_cubic_const = r'''
\s+(?P<index1>\d+)
\s+(?P<index2>\d+)
\s+(?P<index3>\d+)
\s+(?P<const1>[+-]?\d+\.\d+)
\s+(?P<const2>[+-]?\d+\.\d+)
\s+(?P<const3>[+-]?\d+\.\d+)
'''  
pattern_quartic_const = r'''
\s+(?P<index1>\d+)
\s+(?P<index2>\d+)
\s+(?P<index3>\d+)
\s+(?P<index4>\d+)
\s+(?P<const1>[+-]?\d+\.\d+)
\s+(?P<const2>[+-]?\d+\.\d+)
\s+(?P<const3>[+-]?\d+\.\d+)
'''  
pattern_Coriolis_const = r'''
\s+(?P<axes>\w)
\s+(?P<index1>\d+)
\s+(?P<index2>\d+)
\s+(?P<const>[+-]?\d+\.\d+)
'''  
pattern_Rotational_const = r'''
Rotat.Constants:
\s+(?P<const1>[+-]?\d+\.\d+)
\s+(?P<const2>[+-]?\d+\.\d+)
\s+(?P<const3>[+-]?\d+\.\d+)
'''  
#
# DEFINITION OF LOCAL FUNCTIONS
# -----------------------------
def multinomial(list_int):
    '''
    For a given list of integer, computes the associated multinomial
    '''
    tmp      = list(set(list_int))
    multinom = 1
    for i in range(len(tmp)):
        j = list_int.count(tmp[i])
        multinom = multinom*math.factorial(j)
    return multinom
#
# BUILDS THE PARSER
# -----------------
parser = argparse.ArgumentParser(description=DESCR)
parser.add_argument("-log",action='store',type=argparse.FileType(mode="r"),help="Gaussian log file",required=True)
parser.add_argument("-out",action='store',type=argparse.FileType(mode="w"),help="FCIDmp file",required=True)
parser.add_argument("-nvib",action='store',type=int,help="Number of normal modes",required=True)
parser.add_argument("-nofact",action='store_true',help="Don't normalize by the multinomial factor")
parser.add_argument("-DoCor",action='store_true',help="Extract also Coriolis coupling from GAUSSIAN .log file")
options    = parser.parse_args()
log_anharm = options.log
out_file   = options.out
nvib       = options.nvib
DoFact     = options.nofact
DoCor      = options.DoCor
#
# PRINTS THE HEADER
# -----------------
head_template = ''' &FCI NORB= {nvib:d},NELEC=-6,MS2= 2,
  ORBSYM={orbsym:s}
  ISYM=1
 &END
'''
orbsym = '1,'*(nvib-1) + '1'
header = head_template.format(nvib=nvib,orbsym=orbsym)
out_file.write(header)
#
# LOOK FOR ANHARMONIC FORCE FIELD
# -------------------------------
FndCub  = False
FndQrt  = False
FndHrm  = False
FndCor  = False
Lst_Hrm = []
Lst_Cub = []
Lst_Qrt = []
Lst_Rot = []
for i in log_anharm:
    match_Hrm = re.search(pattern_harmonic,i)
    match_Cub = re.search(pattern_cubic,i)
    match_Qrt = re.search(pattern_quartic,i)
    match_Cor = re.search(pattern_Coriolis,i)
    # Looks for the relevant section
    # This part is highly dependent on the version of GAUSSIAN. Might be generalized
    if (match_Cor):
        FndCor = True
        XCor = []
        YCor = []
        ZCor = []
    if (match_Hrm):
        FndHrm = True
    if (match_Cub):
        FndCor = False
        FndCub = True
    if (match_Qrt):
        FndCub = False
        FndQrt = True
    # Looks for the specific patterns
    LJunk0 = re.search(pattern_harmonic_const,i,re.DOTALL|re.VERBOSE)
    LJunk1 = re.search(pattern_cubic_const,i,re.DOTALL|re.VERBOSE)
    LJunk2 = re.search(pattern_quartic_const,i,re.DOTALL|re.VERBOSE)
    LJunk3 = re.search(pattern_Coriolis_const,i,re.DOTALL|re.VERBOSE)
    LJunk4 = re.search(pattern_Rotational_const,i,re.DOTALL|re.VERBOSE)
    if (FndHrm and LJunk0):
        Lst_Hrm.append([int(LJunk0.group('index')),float(LJunk0.group('harm'))])
    if (FndCub and LJunk1):
        Indx1   = int(LJunk1.group('index1'))
        Indx2   = int(LJunk1.group('index2'))
        Indx3   = int(LJunk1.group('index3'))
        LstJnk  = [Indx1,Indx2,Indx3]
        LstJnk2 = [LstJnk,float(LJunk1.group('const1'))]
        Lst_Cub.append(LstJnk2)
    if (FndQrt and LJunk2):
        Indx1   = int(LJunk2.group('index1'))
        Indx2   = int(LJunk2.group('index2'))
        Indx3   = int(LJunk2.group('index3'))
        Indx4   = int(LJunk2.group('index4'))
        LstJnk  = [Indx1,Indx2,Indx3,Indx4]
        LstJnk2 = [LstJnk,float(LJunk1.group('const1'))]
        Lst_Qrt.append(LstJnk2)
    if (FndCor and LJunk3 and DoCor):
        Letter  = LJunk3.group('axes')
        Indx1   = int(LJunk3.group('index1'))
        Indx2   = int(LJunk3.group('index2'))
        LstJnk  = [Indx1,Indx2]
        if(Letter == 'x'):
            XCor.append([LstJnk,float(LJunk3.group('const'))])
        elif(Letter == 'y'):
            YCor.append([LstJnk,float(LJunk3.group('const'))])
        elif(Letter == 'z'):
            ZCor.append([LstJnk,float(LJunk3.group('const'))])
        else:
            print('Unrecognized symbol for axes in zeta Coriolis matrix')
            sys.exit()
    elif (LJunk4 and len(Lst_Rot) == 0 and DoCor):
        Junk1 = float(LJunk4.group('const1'))
        Junk2 = float(LJunk4.group('const2'))
        Junk3 = float(LJunk4.group('const3'))
        Lst_Rot.append(Junk1)
        Lst_Rot.append(Junk2)
        Lst_Rot.append(Junk3)
# Post-processing of Coriolis couplings
if(DoCor):
    Zeta_mat = np.zeros((3,nvib,nvib))
    ICont    = 0
    for Lst in (XCor,YCor,ZCor):
        for j in range(len(Lst)):
            Indx1 = Lst[j][0][0] - 1
            Indx2 = Lst[j][0][1] - 1
    	    # NOTE : the Coriolis Zeta matrix is antysimmetric
            #Zeta_mat[Indx1,Indx2] = Zeta_mat[Indx1,Indx2] + Lst_Rot[ICont]*Lst[j][1]
            #Zeta_mat[Indx2,Indx1] = Zeta_mat[Indx2,Indx1] - Lst_Rot[ICont]*Lst[j][1]
            Zeta_mat[ICont,Indx1,Indx2] =  Lst[j][1]
            Zeta_mat[ICont,Indx2,Indx1] = -Lst[j][1]
        ICont = ICont + 1
#
# PRINTS OUT THE RESULTS
# ----------------------
# NOTE: the factorial factor is added here!
str_qdr = '  {value:19.12E}  {index1:d}  {index2:d}  0  0  0  0 \n'
str_cub = '  {value:19.12E}  {index1:d}  {index2:d}  {index3:d}  0  0  0  \n'
str_qrt = '  {value:19.12E}  {index1:d}  {index2:d}  {index3:d}  {index4:d}  0  0 \n'
for i in range(len(Lst_Hrm)):
    out_file.write(str_qdr.format(value=Lst_Hrm[i][1]/4.,index1=Lst_Hrm[i][0],index2=Lst_Hrm[i][0]))
    out_file.write(str_qdr.format(value=-Lst_Hrm[i][1]/4.,index1=-Lst_Hrm[i][0],index2=-Lst_Hrm[i][0]))
for i in range(len(Lst_Cub)):
    if(not DoFact):
        fact = multinomial(Lst_Cub[i][0])*2.*math.sqrt(2.)
    else:
        fact = 2.*math.sqrt(2.)
    out_file.write(str_cub.format(value=Lst_Cub[i][1]/fact,index1=Lst_Cub[i][0][0],index2=Lst_Cub[i][0][1],index3=Lst_Cub[i][0][2]))
for i in range(len(Lst_Qrt)):
    if(not DoFact):
        fact = multinomial(Lst_Qrt[i][0])
    else:
        fact = 1.
    fact *= 4.
    out_file.write(str_qrt.format(value=Lst_Qrt[i][1]/fact,index1=Lst_Qrt[i][0][0],index2=Lst_Qrt[i][0][1],index3=Lst_Qrt[i][0][2],index4=Lst_Qrt[i][0][3]))
if(DoCor):
    for i in range(nvib):
        for j in range(nvib):
            for k in range(nvib):
                for l in range(nvib):
                    #if (abs(Zeta_mat[i,j]) > 1.0E-10 and abs(Zeta_mat[k,l]) > 1.0E-10):
                    str_Cor = '  {value:19.12E}  {index1:d}  {index2:d}  {index3:d}  {index4:d} \n'
                    Junk  = Zeta_mat[0,i,j]*Zeta_mat[0,k,l]*Lst_Rot[0] + Zeta_mat[1,i,j]*Zeta_mat[1,k,l]*Lst_Rot[1] + Zeta_mat[2,i,j]*Zeta_mat[2,k,l]*Lst_Rot[2]
                    XJunk = Junk*math.sqrt((Lst_Hrm[j][1]*Lst_Hrm[l][1])/(Lst_Hrm[i][1]*Lst_Hrm[k][1]))
                    #XJunk = Zeta_mat[i,j]*Zeta_mat[k,l]
                    out_file.write(str_Cor.format(value=XJunk/4.0,index1=i+1,index2=-j-1,index3=k+1,index4=-l-1))
out_file.close()

