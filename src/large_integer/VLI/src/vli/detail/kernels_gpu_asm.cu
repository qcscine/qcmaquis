



#include "vli/detail/kernels_gpu_asm.h" 


void add192_192(unsigned long int* x /* shared */, unsigned long int const* y /* global */){
   //asm("add.u64 %0, %0, %1  ; \n\t" :"+l"(*x):"l"(*y));
   unsigned int  i,j,k;
   asm("add.u32 %0, %1, %2;" : "=r"(i) : "r"(j), "r"(k));

}