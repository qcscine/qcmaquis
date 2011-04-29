
# Makefile

all :
	make clean gpu 

clean:
	rm $(PWD)/OBJECT_gpu/*; rm -f exe_*  

include Make-config 

D=.
S=.
O=$S/OBJECT_$N
NOM=exe
N=gpu
PROGRAM=$D

OBJ3    = $O/kernels_gpu.o $O/main.o

run :
	$(PROGRAM)/$(NOM)_$(N)

targets_gpu : $(OBJ3)
	$(MPICC) $(OBJ3) $(LINK)  $(NOM) -L$(PATH_LIB_CUBLAS)$(LIB) $(PG)
	mv $(NOM) $(PROGRAM)/$(NOM)_$N



#=============================================
#
#Target !
#=============================================


$O/main.o : 
	$(MPICC) -c -D$(PREPROC) $(OPTIONS) $(PATH_INC) 	main.cpp
	mv main.o $O

$O/kernels_gpu.o : 
	$(NVCC) -c $(OPTIONS_NVCC) $(PATH_INC) 	detail/kernels_gpu.cu
	mv kernels_gpu.o $O
