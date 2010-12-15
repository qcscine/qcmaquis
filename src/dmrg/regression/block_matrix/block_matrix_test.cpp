#include <iostream>
#include <vector>

using namespace std;

#include "dense_matrix/vector_interface.hpp"

#include "dense_matrix/dense_matrix.hpp"
#include "dense_matrix/matrix_interface.hpp"
#include "dense_matrix/resizable_matrix_interface.hpp"
#include "dense_matrix/matrix_algorithms.hpp"
#include "dense_matrix/diagonal_matrix.h"
typedef blas::dense_matrix<double> Matrix;
typedef blas::associated_diagonal_matrix<Matrix>::type DiagMatrix;

#include "block_matrix/block_matrix.h"
#include "block_matrix/block_matrix_algorithms.h"

int main()
{
    {
        typedef U1 grp;
        
        Index<grp> physical;
        physical.push_back(std::make_pair(1, 1));
        physical.push_back(std::make_pair(-2, 1));
        physical.push_back(std::make_pair(-1, 1));
        
        cout << physical << endl;
        physical.sort();
        cout << physical << endl;
        
        cout << physical.position(0) << endl;
        
        physical.insert(make_pair(0,3));
        cout << physical << endl;
        
        cout << physical.position(-1) << endl;
    }
    
    {
        typedef U1 grp;
        
        Matrix foo(10,10);
        block_matrix<Matrix, grp> test;
        test.insert_block(boost::tuples::make_tuple(foo, 0, 0));
        
        cout << test.description() << endl;
    }
    
    if (true)
    {
        typedef U1 grp;
        
        std::vector<grp::charge> c;
        c.push_back(-1); c.push_back(0); c.push_back(1);
        
        std::vector<std::size_t> s(3, 2);
        Index<grp> rows(c, s);
        
        c.clear();
        c.push_back(1); c.push_back(0); c.push_back(-1);
        s.clear();
        s.push_back(2); s.push_back(3); s.push_back(4);
        Index<grp> cols(c, s);
        grp::get_map(cols.charges());
        
        block_matrix<Matrix, grp> m1(rows, cols), m2(cols, rows), m3;
        
        cout << m1.description() << endl << m2.description() << endl;
        cout << m1.description() << endl << m2.description() << endl;
        gemm(m1, m2, m3);
        
        for (std::size_t k = 0; k < m3.n_blocks(); ++k)
            for (std::size_t l = 0; l < num_rows(m3[k]); ++l)
                for (std::size_t r = 0; r < num_columns(m3[k]); ++r)
                    m3[k](l,r) = drand48();
        
        cout << m3 << endl;
        
        block_matrix<Matrix, grp> U, V;
        block_matrix<DiagMatrix, grp> S;
        svd(m3, U, V, S);
        
        //cout << S << endl;
        //cout << U << endl;
        
        block_matrix<Matrix, grp> t;
        gemm(U, S, t);
        gemm(t, V, U);
        
        cout << U << endl;
        
        qr(m3, U, V);
    }
}
