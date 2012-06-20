
#include <cmath>
#include <iterator>
#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

#include "alps/numeric/matrix/matrix.hpp"
#include "alps/numeric/matrix/matrix_interface.hpp"
#include "alps/numeric/matrix/resizable_matrix_interface.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "alps/numeric/matrix/matrix_algorithms.hpp"
#include "alps/numeric/matrix/matrix_blas.hpp"
#include "alps/numeric/diagonal_matrix.hpp"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mpstensor.h"


typedef U1 SymmGroup;
typedef alps::numeric::matrix<double> matrix;


int main() {
    
    int Nrep = 10;
    int M = 50;
    
    // Bosons with Nmax=2
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(2, 1));
    
    for (int i=0; i<Nrep; ++i) {
        // create random left_i
        Index<SymmGroup> left_i;
        size_t max_left = long(dmrg_random::uniform()*10)+1;
        for (size_t k = 0; k<max_left; ++k)
            left_i.insert(std::make_pair(k, long(dmrg_random::uniform()*M)+1));
        
        // create consistent random right_i
        Index<SymmGroup> right_i = phys * left_i;
        for (long k=right_i.size()-1; k>=0; --k) {
            if (dmrg_random::uniform() < 0.2)
                right_i.erase(right_i.begin()+k);
            else
                right_i[k].second = long(dmrg_random::uniform()*M)+1);
        }
        
        MPSTensor<matrix, SymmGroup> m1(phys, left_i, right_i), m2(phys, left_i, right_i);
        
        double norm = m1.scalar_norm();
        double overlap = m1.scalar_overlap(m2);
        
    }

}
