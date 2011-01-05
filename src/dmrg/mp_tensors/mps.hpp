#include "mp_tensors/mps.h"

#include "contractions.h"

#include <boost/math/special_functions/binomial.hpp>

template<class Matrix>
void mps_init(MPS<Matrix, NullGroup> & mps, size_t Mmax, Index<NullGroup> const & phys)
{
    std::size_t L = mps.length();
    std::vector<std::size_t> bond_sizes(L+1, 1);
    for (std::size_t k = 1; k < L+1; ++k)
        bond_sizes[k] = std::min(Mmax, 2*bond_sizes[k-1]);
    bond_sizes[L] = 1;
    for (std::size_t k = L; k > 0; --k) {
        bond_sizes[k-1] = std::min(bond_sizes[k-1], 2*bond_sizes[k]);
    }
    std::copy(bond_sizes.begin(), bond_sizes.end(), std::ostream_iterator<int>(cout, " "));
    cout << endl;
    
    for (int i = 0; i < L; ++i)
    {
        Index<NullGroup> li; li.insert(std::make_pair(NullGroup::Plus, bond_sizes[i]));
        Index<NullGroup> ri; ri.insert(std::make_pair(NullGroup::Plus, bond_sizes[i+1]));
        
        mps[i] = MPSTensor<Matrix, NullGroup>(phys, li, ri);
    }
}

template<class T>
T tri_min(T a, T b, T c)
{
    return std::min(std::min(a, b),
                    std::min(a, c));
}

template<class Matrix>
void mps_init(MPS<Matrix, U1> & mps, size_t Mmax, Index<U1> const & phys)
{
    std::size_t L = mps.length();
    
    Index<U1> triv;
    triv.insert( std::make_pair(U1::SingletCharge, 1) );
    
    std::vector<Index<U1> > left_allowed(L+1), right_allowed(L+1), allowed(L+1);
    left_allowed[0] = triv;
    right_allowed[L] = triv;
    
    for (int i = 1; i < L+1; ++i)
        left_allowed[i] = phys * left_allowed[i-1];
    for (int i = L-1; i >= 0; --i)
        right_allowed[i] = phys * right_allowed[i+1];
    
    for (int i = 0; i < L+1; ++i) {
        allowed[i] = common_subset(left_allowed[i], right_allowed[i]);
        for (typename Index<U1>::iterator it = allowed[i].begin();
             it != allowed[i].end(); ++it)
            it->second = tri_min(Mmax,
                                 left_allowed[i].size_of_block(it->first),
                                 right_allowed[i].size_of_block(it->first));
    }
    
    for (int i = 0; i < L; ++i)
        mps[i] = MPSTensor<Matrix, U1>(phys, allowed[i], allowed[i+1]);
    cout << mps.description() << endl;
}

template<class Matrix, class SymmGroup>
std::string MPS<Matrix, SymmGroup>::description() const
{
    std::ostringstream oss;
    for (int i = 0; i < length(); ++i)
    {
        oss << "MPS site " << i << endl;
        oss << (*this)[i].row_dim() << endl;
        oss << "Sum: " << (*this)[i].row_dim().sum_of_sizes() << endl;
        oss << (*this)[i].col_dim() << endl;
        oss << "Sum: " << (*this)[i].col_dim().sum_of_sizes() << endl;
    }
    return oss.str();
}

template<class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup>::MPS(size_t L, size_t Mmax, Index<SymmGroup> phys)
: std::vector<MPSTensor<Matrix, SymmGroup> >(L)
{
    mps_init<Matrix>(*this, Mmax, phys);
    
    for (int i = 0; i < L; ++i)
        (*this)[i].normalize_left(SVD);
    this->canonize_left();
}

template<class Matrix, class SymmGroup>
Boundary<Matrix, SymmGroup>
MPS<Matrix, SymmGroup>::start_mtx() const
{
    Index<SymmGroup> i; i.insert(std::make_pair(NullGroup::Plus, 1));
    Boundary<Matrix, SymmGroup> ret(i, i, 1);
    ret(0, std::make_pair(SymmGroup::SingletCharge, 0), std::make_pair(SymmGroup::SingletCharge, 0)) = 1;
    return ret;
}

template<class Matrix, class SymmGroup>
typename Matrix::value_type
MPS<Matrix, SymmGroup>::canonize_left()
{
    block_matrix<Matrix, SymmGroup> t;
    for (int i = 0; i < length(); ++i) {
        t = (*this)[i].normalize_left(SVD);
        if (i < length()-1)
            (*this)[i+1].multiply_from_left(t);
    }
    return trace(t);
}

template<class Matrix, class SymmGroup>
typename Matrix::value_type
MPS<Matrix, SymmGroup>::canonize_right()
{
    block_matrix<Matrix, SymmGroup> t;
    for (int i = length()-1; i >= 0; --i) {
        t = (*this)[i].normalize_right(SVD);
        if (i > 0)
            (*this)[i-1].multiply_from_right(t);
    }
    return trace(t);
}

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::normalize_left()
{
    canonize_left();
}

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::normalize_right()
{
    canonize_right();
}

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::canonize(std::size_t center)
{
    for (int i = 0; i < center; ++i)
    {
        block_matrix<Matrix, SymmGroup> t = (*this)[i].normalize_left(SVD);
        (*this)[i+1].multiply_from_left(t);
    }
    
    for (int i = length()-1; i > center; --i)
    {
        block_matrix<Matrix, SymmGroup> t = (*this)[i].normalize_right(SVD);
        (*this)[i-1].multiply_from_right(t);
    }
}

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::grow_l2r_sweep(MPOTensor<Matrix, SymmGroup> const & mpo,
                                            Boundary<Matrix, SymmGroup> const & left,
                                            Boundary<Matrix, SymmGroup> const & right,
                                            std::size_t l, double alpha,
                                            double cutoff,
                                            std::size_t Mmax)
{
    MPSTensor<Matrix, SymmGroup> new_mps = contraction::predict_new_state_l2r_sweep((*this)[l], mpo, left, right, alpha, cutoff, Mmax);
    
    (*this)[l+1] = contraction::predict_lanczos_l2r_sweep((*this)[l+1],
                                                          (*this)[l], new_mps);
    (*this)[l] = new_mps;
}
