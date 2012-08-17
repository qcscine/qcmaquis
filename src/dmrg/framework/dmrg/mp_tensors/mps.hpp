/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/utils/logger.h"
#include "dmrg/mp_tensors/mps.h"
#include "contractions.h"
#include <boost/math/special_functions/binomial.hpp>

#include <limits>

template<class Matrix, class SymmGroup>
std::string MPS<Matrix, SymmGroup>::description() const
{
    std::ostringstream oss;
    for (int i = 0; i < length(); ++i)
    {
        oss << "MPS site " << i << std::endl;
        oss << (*this)[i].row_dim() << std::endl;
        oss << "Sum: " << (*this)[i].row_dim().sum_of_sizes() << std::endl;
        oss << (*this)[i].col_dim() << std::endl;
        oss << "Sum: " << (*this)[i].col_dim().sum_of_sizes() << std::endl;
    }
    return oss.str();
}

template<class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup>::MPS()
: canonized_i(std::numeric_limits<size_t>::max())
{ }

template<class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup>::MPS(size_t L)
: canonized_i(std::numeric_limits<size_t>::max())
, data_(L)
{ }

template<class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup>::MPS(size_t L,
                            size_t Mmax,
                            Index<SymmGroup> phys,
                            typename SymmGroup::charge right_end,
                            mps_initializer<Matrix, SymmGroup> & init)
: canonized_i(std::numeric_limits<size_t>::max())
, data_(L)
{
    init(*this, Mmax, phys, right_end);
    
    // MD: this is actually important
    for (int i = 0; i < L; ++i)
        (*this)[i].normalize_left(DefaultSolver());

    this->normalize_left();
}

template<class Matrix, class SymmGroup>
typename MPS<Matrix, SymmGroup>::value_type const & MPS<Matrix, SymmGroup>::operator[](size_t i) const
{ return data_[i]; }

template<class Matrix, class SymmGroup>
typename MPS<Matrix, SymmGroup>::value_type& MPS<Matrix, SymmGroup>::operator[](size_t i)
{
    if (i != canonized_i)
        canonized_i=std::numeric_limits<size_t>::max();
    return data_[i];
}

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::resize(size_t L)
{
    // if canonized_i < L and L < current L, we could conserve canonized_i
    canonized_i=std::numeric_limits<size_t>::max();
    data_.resize(L);
}

template<class Matrix, class SymmGroup>
size_t MPS<Matrix, SymmGroup>::canonization(bool search) const
{
    if (!search)
        return canonized_i;
    
    size_t center = ((*this)[0].isleftnormalized()) ? 1 : 0;
    for (size_t i=1; i<length(); ++i) {
        if (!(*this)[i].isnormalized() && center != i) {
            canonized_i = std::numeric_limits<size_t>::max();
            return canonized_i;
        } else if ((*this)[i].isleftnormalized() && center == i)
            center = i+1;
        else if ((*this)[i].isleftnormalized()) {
            canonized_i = std::numeric_limits<size_t>::max();
            return canonized_i;
        }
    }
    if (center == length())
        center = length()-1;
    
    canonized_i = center;
    return canonized_i;
}

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::normalize_left()
{
    canonize(length()-1);
    // now state is: A A A A A A M
    block_matrix<Matrix, SymmGroup> t = (*this)[length()-1].normalize_left(DefaultSolver());
    // now state is: A A A A A A A
    canonized_i = length()-1;
}

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::normalize_right()
{
    canonize(0);
    // now state is: M B B B B B B
    block_matrix<Matrix, SymmGroup> t = (*this)[0].normalize_right(DefaultSolver());
    // now state is: B B B B B B B
    canonized_i = 0;
}

// input:  M  M  M  M  M  M  M
//  (idx)        c
// output: A  A  M  B  B  B  B
template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::canonize(std::size_t center, DecompMethod method)
{
    if (canonized_i == center)
        return;
    
    if (canonized_i < center)
        move_normalization_l2r(canonized_i, center, method);
    else if (canonized_i < length())
        move_normalization_r2l(canonized_i, center, method);
    else {
        move_normalization_l2r(0, center, method);
        move_normalization_r2l(length()-1, center, method);
    }
    canonized_i = center;
}

// input:  M  M  M  M  M  M  M
//  (idx)     p1       p2
// output: M  A  A  A  M  M  M
template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::move_normalization_l2r(size_t p1, size_t p2, DecompMethod method)
{
    size_t tmp_i = canonized_i;
    for (int i = p1; i < std::min(p2, length()); ++i)
    {
        if ((*this)[i].isleftnormalized())
            continue;
        block_matrix<Matrix, SymmGroup> t = (*this)[i].normalize_left(method);
        if (i < length()-1) {
            (*this)[i+1].multiply_from_left(t);
            (*this)[i+1].divide_by_scalar((*this)[i+1].scalar_norm());
        }
    }
    if (tmp_i == p1)
        canonized_i = p2;
    else
        canonized_i = std::numeric_limits<size_t>::max();
}

// input:  M  M  M  M  M  M  M
//  (idx)     p2       p1
// output: M  M  B  B  B  M  M
template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::move_normalization_r2l(size_t p1, size_t p2, DecompMethod method)
{
    size_t tmp_i = canonized_i;
    for (int i = p1; i > static_cast<int>(std::max(p2, size_t(0))); --i)
    {
        if ((*this)[i].isrightnormalized())
            continue;
        block_matrix<Matrix, SymmGroup> t = (*this)[i].normalize_right(method);
        if (i > 0) {
            (*this)[i-1].multiply_from_right(t);
            (*this)[i-1].divide_by_scalar((*this)[i-1].scalar_norm());
        }
    }
    if (tmp_i == p1)
        canonized_i = p2;
    else
        canonized_i = std::numeric_limits<size_t>::max();
}

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::grow_l2r_sweep(MPOTensor<Matrix, SymmGroup> const & mpo,
                                            Boundary<Matrix, SymmGroup> const & left,
                                            Boundary<Matrix, SymmGroup> const & right,
                                            std::size_t l, double alpha,
                                            double cutoff, std::size_t Mmax,
                                            Logger & logger)
{ // canonized_i invalided through (*this)[] 
    MPSTensor<Matrix, SymmGroup> new_mps =
    contraction::predict_new_state_l2r_sweep((*this)[l], mpo, left, right, alpha, cutoff, Mmax, logger);
    (*this)[l+1] = contraction::predict_lanczos_l2r_sweep((*this)[l+1],
                                                          (*this)[l], new_mps);
    (*this)[l] = new_mps;
}

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::grow_r2l_sweep(MPOTensor<Matrix, SymmGroup> const & mpo,
                                            Boundary<Matrix, SymmGroup> const & left,
                                            Boundary<Matrix, SymmGroup> const & right,
                                            std::size_t l, double alpha,
                                            double cutoff, std::size_t Mmax,
                                            Logger & logger)
{ // canonized_i invalided through (*this)[] 
    MPSTensor<Matrix, SymmGroup> new_mps =
    contraction::predict_new_state_r2l_sweep((*this)[l], mpo, left, right, alpha, cutoff, Mmax, logger);
    
    (*this)[l-1] = contraction::predict_lanczos_r2l_sweep((*this)[l-1],
                                                          (*this)[l], new_mps);
    (*this)[l] = new_mps;
}

template<class Matrix, class SymmGroup>
Boundary<Matrix, SymmGroup>
MPS<Matrix, SymmGroup>::left_boundary() const
{
    Index<SymmGroup> i = (*this)[0].row_dim();
    Boundary<Matrix, SymmGroup> ret(i, i, 1);

    for(std::size_t k(0); k < ret[0].n_blocks(); ++k)
       maquis::dmrg::detail::left_right_boundary_init(ret[0][k]);

    return ret;
}

template<class Matrix, class SymmGroup>
Boundary<Matrix, SymmGroup>
MPS<Matrix, SymmGroup>::right_boundary() const
{
    Index<SymmGroup> i = (*this)[length()-1].col_dim();
    Boundary<Matrix, SymmGroup> ret(i, i, 1);

//    Original
//    for(typename Index<SymmGroup>::basis_iterator it = i.basis_begin(); !it.end(); ++it)
//        ret(0,*it,*it) = 1;

    for(std::size_t k(0); k < ret[0].n_blocks(); ++k)
        maquis::dmrg::detail::left_right_boundary_init(ret[0][k]);

    return ret;
}

#ifdef HAVE_ALPS_HDF5

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::load(alps::hdf5::archive & ar)
{
    canonized_i = std::numeric_limits<size_t>::max();
    ar >> alps::make_pvp("MPS", data_);
}

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::save(alps::hdf5::archive & ar) const
{
    ar << alps::make_pvp("MPS", data_);
}

#endif


template <class Matrix, class SymmGroup>
void check_equal_mps (MPS<Matrix, SymmGroup> const & mps1, MPS<Matrix, SymmGroup> const & mps2)
{
    // Length
    if (mps1.length() != mps2.length())
        throw std::runtime_error("Length doesn't match.");
    
    for (int i=0; i<mps1.length(); ++i)
        try {
            mps1[i].check_equal(mps2[i]);
        } catch (std::exception & e) {
            maquis::cerr << "Problem on site " << i << ":" << e.what() << std::endl;
            exit(1);
        }
}


