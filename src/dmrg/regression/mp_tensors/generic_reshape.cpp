#include <iostream>

#include <vector>
#include <iterator>
#include <algorithm>
#include <cstdlib>

#include "dmrg/block_matrix/detail/alps_matrix.hpp"
typedef alps::numeric::matrix<double> Matrix;


#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/multi_index.h"
#include "dmrg/mp_tensors/reshapes.h"
#include "dmrg/mp_tensors/generic_reshape.h"

typedef U1 symm;

std::ostream& operator<< (std::ostream& os, std::pair<symm::charge, std::size_t> const& p)
{
    os << "(" << p.first << " : " << p.second << ")";
    return os;
}

std::ostream& operator<< (std::ostream& os, index_product_iterator<symm>::value_type const& v)
{
    //std::copy(v.begin(), v.end(), std::ostream_iterator<std::pair<symm::charge, std::size_t> >(os, " "));
    for (int i=0; i<v.size(); ++i)
        os << v[i] << " ";
    return os;
}

std::ostream& operator<< (std::ostream& os, std::pair<MultiIndex<symm>::coord_t, MultiIndex<symm>::coord_t> const& p)
{
    os << p.first << ", " << p.second;
    return os;
}


void test_reshape() {
    
    Index<symm> phys_i, alpha_i, beta_i;
    phys_i.insert(std::make_pair(0, 1));
    phys_i.insert(std::make_pair(1, 1));
    
    // create alpha, beta indexes
    for (size_t n=0; n<100; ++n)
        alpha_i.insert(std::make_pair(n, 4/(n+1) + 1));

    for (size_t n=0; n<400; ++n)
        beta_i.insert(std::make_pair(n, 3/(n+1) + 1));

//    beta_i = phys_i * alpha_i;
    
    MultiIndex<symm> midx;
    MultiIndex<symm>::index_id alpha = midx.insert_index(alpha_i);
    MultiIndex<symm>::index_id beta = midx.insert_index(beta_i);
    MultiIndex<symm>::index_id sigma = midx.insert_index(phys_i);
    
    
    std::vector< std::pair<MultiIndex<symm>::index_id, bool> > left_paired_rows(2), left_paired_cols(1);
    left_paired_rows[0] = std::make_pair(sigma, true); left_paired_rows[1] = std::make_pair(alpha, true);
    left_paired_cols[0] = std::make_pair(beta, true);
    MultiIndex<symm>::set_id left_paired = midx.create_set(left_paired_rows, left_paired_cols);
    
    std::vector< std::pair<MultiIndex<symm>::index_id, bool> > right_paired_rows(1), right_paired_cols(2);
    right_paired_rows[0] = std::make_pair(alpha, true);
    right_paired_cols[0] = std::make_pair(sigma, false); right_paired_cols[1] = std::make_pair(beta, true);
    MultiIndex<symm>::set_id right_paired = midx.create_set(right_paired_rows, right_paired_cols);
    
    
    // fill block_matrix
    maquis::cout << "* filling matrices" << std::endl;
    Index<symm> left_i = phys_i*alpha_i;
    Index<symm> right_i = adjoin(phys_i)*beta_i;
    
    maquis::cout << " - allocating matrix" << std::endl;
    block_matrix<Matrix, symm> m1; // left-paired
    size_t max_i = std::min(left_i.size(), beta_i.size());
    for(size_t i=0; i<max_i; ++i)
        m1.insert_block(Matrix(left_i[i].second, beta_i[i].second), left_i[i].first, beta_i[i].first);
    
    
    maquis::cout << " - filling values" << std::endl;
    m1.generate(drand48);

//    maquis::cout << "alpha basis: " << alpha_i << std::endl;
//    maquis::cout << "beta basis: " << beta_i << std::endl;
//    maquis::cout << "sigma basis: " << phys_i << std::endl;
//    maquis::cout << "left_paired basis: " << left_i << std::endl;
//    maquis::cout << "right_paired basis: " << right_i << std::endl;
//    
//    maquis::cout << "left basis: " << m1.left_basis() << std::endl;
//    maquis::cout << "right basis: " << m1.right_basis() << std::endl;
//    
//    for(index_product_iterator<symm> it(midx.begin());
//        it != midx.end();
//        it++)
//    {
//        maquis::cout << *it << " = " << midx.get_coords(left_paired, *it)
//        << " or " << midx.get_coords(right_paired, *it)
//        << std::endl;
//    }

    
    
    maquis::cout << " - copying m1" << std::endl;
    block_matrix<Matrix, symm> m2(m1); // left-paired
    size_t n_blocks = m2.n_blocks();
    while (m2.n_blocks() > n_blocks/2) {
        m2.remove_block( static_cast<size_t>( drand48()*m2.n_blocks() ) );
    }
    
    
    // run reshape and reshape_left_to_right
    block_matrix<Matrix, symm> res1, res1_gen, res2, res2_gen;
    
    maquis::cout << "* standard reshape" << std::endl;
    reshape_left_to_right(phys_i, alpha_i, beta_i, m1, res1);
    reshape_left_to_right(phys_i, alpha_i, beta_i, m2, res2);

    maquis::cout << "* generic reshape" << std::endl;
    reshape(midx, left_paired, right_paired, m1, res1_gen);
    reshape(midx, left_paired, right_paired, m2, res2_gen);

    
    // compare results
    maquis::cout << "* comparing" << std::endl;
    assert( res1.n_blocks() == res1_gen.n_blocks() );
    assert( res1.left_basis() == res1_gen.left_basis() );
    assert( res1.right_basis() == res1_gen.right_basis() );
    for (size_t n=0; n<res1.n_blocks(); ++n)
        assert( equal(res1[n].elements().first, res1[n].elements().first,
                      res1_gen[n].elements().first) );
    
    assert( res2.n_blocks() == res2_gen.n_blocks() );
    assert( res2.left_basis() == res2_gen.left_basis() );
    assert( res2.right_basis() == res2_gen.right_basis() );
    for (size_t n=0; n<res2.n_blocks(); ++n)
        assert( equal(res2[n].elements().first, res2[n].elements().first,
                      res2_gen[n].elements().first) );
   
    
    
    // benchmarks
    maquis::cout << "* benchmarks" << std::endl;
    for (int n=0; n<1000; ++n) {
        reshape_left_to_right(phys_i, alpha_i, beta_i, m1, res1);
        reshape(midx, left_paired, right_paired, m1, res1_gen);
    }
    
}


int main () {
    maquis::cout << "*** Generic reshape" << std::endl;
    test_reshape();
    
    
}

