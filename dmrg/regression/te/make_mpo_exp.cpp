#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/aligned_allocator.h"
#include "types/dense_matrix/algorithms.hpp"
typedef maquis::types::dense_matrix<double> Matrix;

#include <alps/hdf5.hpp>

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/mp_tensors/mps_initializers.h"


#include "dmrg/models/hamiltonian.h"

#include "dmrg/mp_tensors/te.h"
#include "te_utils.hpp"


using namespace app;

typedef U1 grp;

typedef std::vector<MPOTensor<Matrix, grp> > mpo_t;
typedef Boundary<Matrix, grp> boundary_t;

std::ostream& operator<< (std::ostream& os, std::pair<grp::charge, std::size_t> const& p)
{
    os << "(" << p.first << " : " << p.second << ")";
    return os;
}

std::ostream& operator<< (std::ostream& os, index_product_iterator<grp>::value_type const& v)
{
    //std::copy(v.begin(), v.end(), std::ostream_iterator<std::pair<grp::charge, std::size_t> >(os, " "));
    for (int i=0; i<v.size(); ++i)
        os << v[i] << " ";
    return os;
}

std::ostream& operator<< (std::ostream& os, std::pair<MultiIndex<grp>::coord_t, MultiIndex<grp>::coord_t> const& p)
{
    os << p.first << ", " << p.second;
    return os;
}

std::ostream& operator<< (std::ostream& os, MPO<Matrix, grp> const& mpo)
{
    for (int p=0; p<mpo.size(); ++p)
        for (size_t r=0; r<mpo[p].row_dim(); ++r)
            for (size_t c=0; c<mpo[p].col_dim(); ++c)
                if (mpo[p].has(r, c))
                    os << "** Position " << p << " [" << r << "," << c << "]:" << std::endl << mpo[p](r,c);
    return os;
}

Hamiltonian<Matrix, U1> create_H (size_t length, bool with_sign=false)
{
    typedef Hamiltonian<Matrix, U1> ham;        
    typedef ham::hamterm_t hamterm_t;        
    typedef ham::op_t op_t;
    
    double t=1;
    
    op_t ident;
    op_t create, destroy, count;
    
    ident.insert_block(Matrix(1, 1, 1), 0, 0);
    ident.insert_block(Matrix(1, 1, 1), 1, 1);
    
    create.insert_block(Matrix(1, 1, 1), 0, 1);
    destroy.insert_block(Matrix(1, 1, 1), 1, 0);
    
    count.insert_block(Matrix(1, 1, 1), 1, 1);
    
    Index<U1> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    
    std::vector<hamterm_t> terms;
    for (int p=0; p<length-1; ++p) {
        {
            hamterm_t term;
            term.with_sign = with_sign;
            term.fill_operator = ident;
            term.operators.push_back( std::make_pair(p, -t*create) );
            term.operators.push_back( std::make_pair(p+1, destroy) );
            terms.push_back(term);
        }
        {
            hamterm_t term;
            term.with_sign = with_sign;
            term.fill_operator = ident;
            term.operators.push_back( std::make_pair(p, -t*destroy) );
            term.operators.push_back( std::make_pair(p+1, create) );
            terms.push_back(term);
        }
    }
    
    return ham(phys, ident, terms);
}


std::vector<MPO<Matrix, grp> >
getU(std::vector<Hamiltonian<Matrix, grp> > const & split_H, size_t length)
{
//    double alpha = -dt;
    double alpha = 1;
        
    std::vector<MPO<Matrix, grp> > expMPO(split_H.size(), MPO<Matrix, grp>(length));
    for (int i=0; i<split_H.size(); ++i) {
        expMPO[i] = make_exp_mpo(length, split_H[i], alpha);

    }
    return expMPO;
}

int main(int argc, char ** argv)
{
    size_t length = 4;
    
    Hamiltonian<Matrix, grp> H = create_H(length);
    Hamiltonian<Matrix, grp> H_signed = create_H(length, true);
//    maquis::cout << "** Hamiltonian **" << std::endl << H;
    
    std::vector<Hamiltonian<Matrix, grp> > split_H = separate_overlaps(H);
    std::vector<Hamiltonian<Matrix, grp> > split_H_signed = separate_overlaps(H_signed);
    maquis::cout << "length split_H = " << split_H.size() << std::endl;
    split_H.erase(split_H.begin() + 1);
    MPO<Matrix, grp> mpo = make_mpo(length, split_H[0]);
    std::vector<MPO<Matrix, grp> > expMPO = getU(split_H, length);
    std::vector<MPO<Matrix, grp> > expMPO_gen = getU(split_H_signed, length);
    
    maquis::cout << "** MPO **" << std::endl << mpo;
    maquis::cout << "** expMPO **" << std::endl << expMPO[0];
    maquis::cout << "** expMPO_gen **" << std::endl << expMPO_gen[0];
}

