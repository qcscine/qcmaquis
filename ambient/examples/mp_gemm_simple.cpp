#include "ambient/ambient.hpp"
#include "ambient/container/numeric/matrix.hpp"
#include "utils/timings.hpp"

#define N AMBIENT_DEFAULT_IB*4 

ambient::scope::const_iterator where(int i, int j){
    static int nq = ambient::scope::size();
    static int np = 1;
    return (ambient::scope::begin() + (i % np) * nq + j % nq);
}

int main(){
    ambient::numeric::tiles< ambient::numeric::matrix<double> > a(N,N), b(N,N), c(N,N);

    for(int i = 0; i < a.mt; i++)
    for(int j = 0; j < a.nt; j++){
        ambient::actor proc(where(i,j));
        ambient::numeric::fill_random(a.tile(i,j));
        ambient::numeric::fill_random(b.tile(i,j));
    }

    ambient::timer t1("gemm"); t1.begin(); std::cout << "gemm (" << N << ")\n";

    for(int k = 0; k < a.nt; k++)
    for(int j = 0; j < c.nt; j++)
    for(int i = 0; i < c.mt; i++){
        ambient::actor proc(where(i,j));
        ambient::numeric::gemm_fma(a.tile(i,k), b.tile(k,j), c.tile(i,j));
    }

    t1.end();
    return 0;
}

