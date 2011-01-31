
# Makefile

default:
	@echo "Missing parameters :"
	@echo "clean, cpu_gcc "


clean:
	rm $$(PWD)/OBJECT_gpu/*; rm -f exe_*  

include Make-config 

D=.
S=.
O=$S/OBJECT_$N
NOM=exe
N=gpu
PROGRAM=$D

OBJ1    = $O/u1_dmrg.o  
OBJ2    = $O/dgemmcpugpu.o
OBJ3    = $O/kernels.o $O/transposecpu_gpu.o

all :
	make clean gpu run 
run :
	$(PROGRAM)/$(NOM)_$(N)

targets : $(OBJ1)
	$(MPICC) $(OBJ1) $(LINK) $(NOM) -L$(PATH_LIB_CUBLAS) -L$(PATH_LIB_LAPACK) -L$(PATH_LIB_ALPS) -L$(PATH_LIB_BOOST) -L$(PATH_LIB_ARPACK) $(LIB) $(PG)
	mv $(NOM) $(PROGRAM)/$(NOM)_$N_sse
	
targets_bench : $(OBJ2)
	$(MPICC) $(OBJ2) $(LINK) $(NOM) -L$(PATH_LIB_CUBLAS) -L$(PATH_LIB_LAPACK) $(LIB) $(PG)
	mv $(NOM) $(PROGRAM)/$(NOM)_$N

targets_gpu : $(OBJ3)
	$(MPICC) $(OBJ3) $(LINK)  $(NOM) -L/usr/local/cuda/lib /usr/local/lib/libreflapack.a -L/usr/local/cula/lib64 /usr/local/lib/librefblas.a -L/opt/alps/lib -L/opt/boost/lib /Users/Tim/Documents/DMRG_XCODE/dmrg/ARPACK/libarpack_MAC.a -lcudart -lcublas -lcuda -lhdf5 -lalps -lboost -lgfortran
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

$O/transposecpu_gpu.o : 
	$(MPICC) -c $(OPTIONS) $(PATH_INC) 	/Users/Tim/Documents/DMRG_XCODE/dmrg/MPS/trunk/src/dmrg/regression/dense_matrix/transposecpu_gpu.cpp
	mv transposecpu_gpu.o $O

$O/kernels.o : 
	$(NVCC) -c $(OPTIONS_NVCC) $(PATH_INC) 	/Users/Tim/Documents/DMRG_XCODE/dmrg/MPS/trunk/src/dmrg/dense_matrix/gpu/kernels.cu
	mv kernels.o $O
