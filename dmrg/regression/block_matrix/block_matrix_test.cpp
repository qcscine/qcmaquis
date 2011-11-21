#include <iostream>
#include <vector>

using namespace std;

#include "utils/zout.hpp"
#include "types/dense_matrix/vector_interface.hpp"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/diagonal_matrix.h"
#include "types/dense_matrix/dense_matrix_algorithms.h"
typedef maquis::types::dense_matrix<double> Matrix;
typedef maquis::types::associated_real_diagonal_matrix<Matrix>::type DiagMatrix;

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

int main()
{
    {
        typedef U1 grp;
        
        Index<grp> physical;
        physical.push_back(std::make_pair(1, 1));
        physical.push_back(std::make_pair(-2, 1));
        physical.push_back(std::make_pair(-1, 1));
        
        zout << physical << endl;
        physical.sort();
        zout << physical << endl;
        
        zout << physical.position(0) << endl;
        
        physical.insert(make_pair(0,3));
        zout << physical << endl;
        
        zout << physical.position(-1) << endl;
    }
    
    {
        typedef U1 grp;
        
        Matrix foo(10,10);
        block_matrix<Matrix, grp> test;
        test.insert_block(foo, 0, 0);
        
        zout << test.description() << endl;
    }
    
    if (true)
    {
        typedef U1 grp;
        
        Index<grp> rc;
        rc.insert(std::make_pair(-1, 2));
        rc.insert(std::make_pair( 0, 2));
        rc.insert(std::make_pair(+1, 2));
        
        block_matrix<Matrix, grp> m1(rc, rc), m2(rc, rc), m3;
        
        zout << m1.description() << endl << m2.description() << endl;
        zout << m1.description() << endl << m2.description() << endl;
        gemm(m1, m2, m3);
        
        for (std::size_t k = 0; k < m3.n_blocks(); ++k)
            for (std::size_t l = 0; l < num_rows(m3[k]); ++l)
                for (std::size_t r = 0; r < num_cols(m3[k]); ++r)
                    m3[k](l,r) = drand48();
        
        zout << m3 << endl;
        
        block_matrix<Matrix, grp> U, V;
        block_matrix<DiagMatrix, grp> S;
        svd(m3, U, V, S);
        
        //zout << S << endl;
        //zout << U << endl;
        
        block_matrix<Matrix, grp> t;
        gemm(U, S, t);
        gemm(t, V, U);
        
        zout << U << endl;
        
        qr(m3, U, V);
    }
}
