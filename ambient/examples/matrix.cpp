#include "ambient/ambient.hpp"
#include "ambient/numeric/matrix.hpp"

int main(){
    using namespace ambient;
    using namespace ambient::numeric;
    typedef tiles<matrix<double> > mtx;

    size_t x = 4096;
    size_t y = 4096;

    mtx pA(x,y);
    mtx pB(x,y);
    mtx pC(x,y);
    mtx pC_orig(x,y);

    generate(pA);
    generate(pB);

    cout << "gemm...\n"; timer t2("gemm"); t2.begin();
    gemm(pA, pB, pC_orig); 
    t2.end();

    cout << "strassen gemm...\n"; timer t1("gemm_strassen"); t1.begin();
    gemm_strassen(std::move(pA), std::move(pB), std::move(pC)); 
    t1.end();

    if(pC == pC_orig) cout << "Success\n";
    return 0;
}

