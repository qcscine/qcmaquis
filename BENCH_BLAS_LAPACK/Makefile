
# Makefile

default:
	@echo "Missing parameters :"
	@echo "clean, ACML, run_ACML ... "


clean:
	rm -f `find . -name \*.o`; rm -f `find . -name \exe_*`  

include Make-config 

D=.
S=.
O=$S/OBJECT_$N
NROWS = 4000 
NCOLS = 2000
PROGRAM=$D
NOM=exe
N=par
N1=GOTO
N2=ACML
N3=INTEL

OBJ1    = $O/main.o $O/CSMatrix.o

lbatch :
	@echo "mprun -np $(NUM)x$(NUMTHREAD) -p $(PART) $(PROGRAM)/$(NOM)_$(N)" | at -m now
all :
	make clean GOTO ACML INTEL
	
run_GOTO :
	for I in 1 2 4 8 16 32 48; do export GOTO_NUM_THREADS=$$I; time taskset -c 1-48 $(PROGRAM)/$(NOM)_$(N1) $(NROWS) $(NCOLS) $$I "GOTO"; done
	
run_ACML :
	for I in 1 2 4 8 16 32 48; do export OMP_NUM_THREADS=$$I; time taskset -c 1-48 $(PROGRAM)/$(NOM)_$(N2) $(NROWS) $(NCOLS) $$I "ACML"; done
	
run_INTEL : 
	export OMP_NUM_THREADS=$(NUMTHREAD); mpiexec -np $(NUM)  $(PROGRAM)/$(NOM)_$(N3) $(NROWS) $(NCOLS) $(NROWSBLOCKS) $(NCOLSBLOCKS)

targets : $(OBJ1)
	$(MPICC) $(OBJ1) $(LINK) $(NOM) -L$(PATH_LIB1)  $(PG)
	mv $(NOM) $(PROGRAM)/$(NOM)_$N

#=============================================
#
#Target !
#=============================================


$O/main.o : 
	$(MPICC) -c $(OPTIONSmpi)	main.cpp
	mv main.o $O

$O/CSMatrix.o : 
	$(MPICC) -c $(OPTIONSmpi)	CSMatrix.cpp
	mv CSMatrix.o $O


