#include "mps.h"

template<class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup>::MPS(size_t L, size_t Mmax, Index<SymmGroup> phys)
: std::vector<MPSTensor<Matrix, SymmGroup> >(L)
{
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
        Index<SymmGroup> li; li.insert(std::make_pair(NullGroup::Plus, bond_sizes[i]));
        Index<SymmGroup> ri; ri.insert(std::make_pair(NullGroup::Plus, bond_sizes[i+1]));
        
        (*this)[i] = MPSTensor<Matrix, SymmGroup>(phys, li, ri);
    }
    
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
