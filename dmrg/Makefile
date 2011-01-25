
# Makefile

default:
	@echo "Missing parameters :"
	@echo "clean, cpu_gcc "


clean:
	rm OBJECT_*/*.o; rm -f exe_*  

include Make-config 

D=.
S=.
O=$S/OBJECT_$N
NOM=exe
N=cpu_gcc
PROGRAM=$D

OBJ1    = $O/u1_dmrg.o  
OBJ2    = $O/dgemmcpugpu.o


all :
	make clean par run 
run :
	$(PROGRAM)/$(NOM)_$(N) params model outfiles

targets : $(OBJ1)
	$(MPICC) $(OBJ1) $(LINK) $(NOM) -L$(PATH_LIB_CUBLAS) -L$(PATH_LIB_LAPACK) -L$(PATH_LIB_ALPS) -L$(PATH_LIB_BOOST) -L$(PATH_LIB_ARPACK) $(LIB) $(PG)
	mv $(NOM) $(PROGRAM)/$(NOM)_$N_sse
	
targets_bench : $(OBJ2)
	$(MPICC) $(OBJ2) $(LINK) $(NOM) -L$(PATH_LIB_CUBLAS) -L$(PATH_LIB_LAPACK) $(LIB) $(PG)
	mv $(NOM) $(PROGRAM)/$(NOM)_$N


#=============================================
#
#Target !
#=============================================


$O/u1_dmrg.o : 
	$(MPICC) -c $(OPTIONS) $(PATH_INC)	u1_dmrg.cpp
	mv u1_dmrg.o $O


$O/dgemmcpugpu.o : 
	$(MPICC) -c $(OPTIONS) $(PATH_INC)	dgemmcpugpu.cpp
	mv dgemmcpugpu.o $O



