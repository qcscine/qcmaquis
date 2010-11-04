#include <iostream>
#include <vector>

using namespace std;

#include <boost/numeric/bindings/std/vector.hpp>
#include "general_matrix.hpp"
typedef blas::general_matrix<double> Matrix;

#include "block_matrix.h"

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
        
        cout << physical.destination(0) << endl;
        
        physical.insert(make_pair(0,3));
        cout << physical << endl;
        
        cout << physical.at(-1) << endl;
    }
    
    {
        typedef U1 grp;
        
        Matrix foo(10,10);
        block_matrix<Matrix, grp> test;
        test.insert_block(boost::tuples::make_tuple(foo, 0, 0));
        
        cout << test.description() << endl;
    }
    
    if (false)
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
        cout << m1.description() << endl << m2.description() << endl;
        gemm(m1, m2, m3);
        
        block_matrix<Matrix, grp> U, S, V;
        svd(m3, U, V, S);
        
        cout << S << endl;
        cout << U << endl;
    }
}
