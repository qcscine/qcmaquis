
# Makefile

default:
	@echo "Missing parameters :"
	@echo "clean, par, par_ACML ... "


clean:
	-rm -f `find . -name \*.o`  

F_LIBS= -ldxml
include Make-config 

D=.
S=.
O=$S/OBJECT_$N
NUM = 16
NUMTHREAD = 1
NROWS = 5000
NCOLS = 4800
NBLOCKS = 32
PROGRAM=$D
NOM=exe
N=par
N2=par_ACML

OBJ1    = $O/main_mpi.o $O/Crand.o $O/CMatrix.o

lbatch :
	@echo "mprun -np $(NUM)x$(NUMTHREAD) -p $(PART) $(PROGRAM)/$(NOM)_$(N)" | at -m now
all :
	make clean par par_ACML
run :
	mpirun -np $(NUM) ./taskit $(PROGRAM)/$(NOM)_$(N) $(NROWS) $(NCOLS) $(NBLOCKS) 
run_ACML :
	export OMP_NUM_THREADS=$(NUMTHREAD); mpirun -np $(NUM) ./taskit $(PROGRAM)/$(NOM)_$(N2) $(NROWS) $(NCOLS) $(NBLOCKS) 

targets : $(OBJ1)
	$(MPICC) $(OBJ1) $(LINK) $(NOM) -L$(PATH_LIB1)  $(PG)
	mv $(NOM) $(PROGRAM)/$(NOM)_$N

#=============================================
#
#Target !
#=============================================


$O/main_mpi.o : 
	$(MPICC) -c $(OPTIONSmpi)	main_mpi.cpp
	mv main_mpi.o $O

$O/Crand.o : 
	$(MPICC) -c $(OPTIONSmpi)	Crand.cpp
	mv Crand.o $O

$O/CMatrix.o : 
	$(MPICC) -c $(OPTIONSmpi)	CMatrix.cpp
	mv CMatrix.o $O

