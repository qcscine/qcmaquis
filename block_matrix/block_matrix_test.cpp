#include <iostream>
#include <vector>

using namespace std;

// #include <alps/numeric/detail/matrix.hpp>
// typedef blas::matrix Matrix;

#include <alps/numeric/detail/general_matrix.hpp>
typedef blas::general_matrix<double> Matrix;

void gemm(Matrix const & A, Matrix const & B, Matrix & C)
{
    A.matrix_right_multiply(B, C, 1);
}

#include "block_matrix.h"

double gemm(double a, double b, double c)
{
    c = a*b;
}

int main()
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
    
    block_matrix<Matrix, grp> m1(rows, cols), m2(transpose(cols), transpose(rows)), m3;
    
    cout << m1.description() << endl << m2.description() << endl;
    match_blocks(m1, m2);
    cout << m1.description() << endl << m2.description() << endl;
    gemm(m1, m2, m3);
}
