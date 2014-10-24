#include "ambient/ambient.hpp"
#include "ambient/container/numeric/matrix.hpp"

int main(){
    using namespace ambient;
    using namespace ambient::numeric;
    typedef tiles<matrix<double> > mtx;

    size_t m = 16384;
    size_t n = 16384;

    mtx pA(m,n);
    mtx pB(n,m);
    mtx pC(m,m);
    mtx pC_orig(m,m);

    generate(pA);
    generate(pB);

    timer t1("tiled gemm"); t1.begin(); cout << "tiled gemm...\n"; 
    gemm(pA, pB, pC); 
    t1.end();

    merge(pA); merge(pB); merge(pC_orig); timer t2("single block gemm"); t2.begin(); cout << "single block gemm...\n"; 
    gemm(pA[0], pB[0], pC_orig[0]);
    ambient::sync(mkl_parallel());
    t2.end(); split(pA); split(pB); split(pC_orig);

    if(pC == pC_orig) cout << "\n";
    return 0;
}

