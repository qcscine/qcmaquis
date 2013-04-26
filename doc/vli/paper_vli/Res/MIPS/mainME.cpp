#include <iostream>
#include <math.h>

double foo( int Order, int var, int mul, int add){
    return 4096*pow((Order+1),2*var)*mul +4095*pow((2*Order+1),var)*add;
}

int main(){

for ( int i=0; i < 15; ++i)
     std::cout << i << " " <<  foo(i,1,32,0) << " " <<  foo(i,2,32,0)<<  " " <<  foo(i,3,32,0) <<  " " <<  foo(i,4,32,0) << std::endl;
    // std::cout << i << " " <<  foo(i,1,1,1) << " " <<  foo(i,2,1,1)<<  " " <<  foo(i,3,1,1) <<  " " <<  foo(i,4,1,1) << std::endl;






}
