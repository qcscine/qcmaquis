#include "ambient/ambient.hpp"
#include "ambient/container/numeric/matrix.hpp"
#include "utils/timings.hpp"

#define GRAN 2
#define N AMBIENT_IB*GRAN*2

namespace ambient { namespace numeric {

    int where(int i, int j){
        static int nq = ambient::scope::size();
        static int np = 1;
        return (i % np) * nq + j % nq;
    }

    template<int G, class Matrix>
    void mp_generate(tiles<Matrix>& a){
        for(int i = 0; i < a.mt/G; i++)
        for(int j = 0; j < a.nt/G; j++){
            ambient::actor proc(ambient::scope::begin() + where(i,j));
            auto tile = a.subset(i*G, j*G, G, G);
            for(int k = 0; k < G*G; k++) fill_random(tile[k]);
        }
    }

    template<int G, class Matrix>
    void mp_gemm(const tiles<Matrix>& a, const tiles<Matrix>& b, tiles<Matrix>& c){
        int mt = c.mt / G;
        int nt = c.nt / G;
        int zt = a.nt / G;
        for(int i = 0; i < mt; i++)
        for(int j = 0; j < nt; j++){
            std::vector<std::pair<tiles<Matrix>*, int> > ctree;
            auto res = c.subset(i*G, j*G, G, G);
            size_t m = res.num_rows();
            size_t n = res.num_cols();
            for(int k = 0; k < zt; k++){
                ambient::actor proc(ambient::scope::begin() + where(k,j));
                ctree.push_back(std::make_pair(new tiles<Matrix>(m,n), where(k,j)));
                gemm(a.subset(i*G, k*G, G, G), 
                     b.subset(k*G, j*G, G, G),
                     *ctree.back().first);
            }
            int root = where(i,j);
            std::sort(ctree.begin(), ctree.end(), [root](const std::pair<tiles<Matrix>*, int>& a, const std::pair<tiles<Matrix>*, int>& b){
                return (ambient::scope::size() + a.second - root) % ambient::scope::size()
                     < (ambient::scope::size() + b.second - root) % ambient::scope::size();
            });
            res = *ambient::reduce(ctree, [](std::pair<tiles<Matrix>*, int>& a, std::pair<tiles<Matrix>*, int>& b){
                ambient::actor proc(ambient::scope::begin() + a.second);
                *a.first += *b.first;
            }).first;
            for(auto x : ctree) delete x.first;
        }
    }
} }

int main(){
    using namespace ambient;
    using namespace ambient::numeric;
    typedef tiles<matrix<double> > mtx;

    mtx pA(N,N); mp_generate<GRAN>(pA);
    mtx pB(N,N); mp_generate<GRAN>(pB);
    mtx pC(N,N);

    ambient::timer t1("gemm"); t1.begin();
    std::cout << "distributed gemm... (dim " << N << ")\n";

    mp_gemm<GRAN>(pA, pB, pC);

    t1.end();
    return 0;
}

