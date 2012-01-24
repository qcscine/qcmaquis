/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#include "utils/zout.hpp"
#include "dmrg/utils/logger.h"

#include "dmrg/mp_tensors/mps.h"

#include "contractions.h"
#include "dmrg/mp_tensors/detail/algorithms_impl.h"
#include <boost/math/special_functions/binomial.hpp>

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
MPS<Matrix, SymmGroup>::MPS(size_t L)
: std::vector<MPSTensor<Matrix, SymmGroup> >(L)
{ }

template<class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup>::MPS(size_t L,
                            size_t Mmax,
                            Index<SymmGroup> phys,
                            typename SymmGroup::charge right_end,
                            mps_initializer<Matrix, SymmGroup> & init)
: std::vector<MPSTensor<Matrix, SymmGroup> >(L)
{
    init(*this, Mmax, phys, right_end);
    
    for (int i = 0; i < L; ++i)
        (*this)[i].normalize_left(SVD);

    this->canonize_left();
}

template<class Matrix, class SymmGroup>
typename Matrix::value_type
MPS<Matrix, SymmGroup>::canonize_left()
{
    block_matrix<Matrix, SymmGroup> t;
    for(int i = 0; i < length(); ++i){
        t = (*this)[i].normalize_left(SVD);
        if(i < length()-1) {
            (*this)[i+1].multiply_from_left(t);
            (*this)[i+1].multiply_by_scalar(1. / (*this)[i+1].scalar_norm());
        }
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
        if (i > 0) {
            (*this)[i-1].multiply_from_right(t);
            (*this)[i-1].multiply_by_scalar(1. / (*this)[i-1].scalar_norm());
        }
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
        (*this)[i+1].multiply_by_scalar(1. / (*this)[i+1].scalar_norm());
    }
    
    for (int i = length()-1; i > center; --i)
    {
        block_matrix<Matrix, SymmGroup> t = (*this)[i].normalize_right(SVD);
        (*this)[i-1].multiply_from_right(t);
        (*this)[i-1].multiply_by_scalar(1. / (*this)[i-1].scalar_norm());
    }
}

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::grow_l2r_sweep(MPOTensor<Matrix, SymmGroup> const & mpo,
                                            Boundary<Matrix, SymmGroup> const & left,
                                            Boundary<Matrix, SymmGroup> const & right,
                                            std::size_t l, double alpha,
                                            double cutoff, std::size_t Mmax,
                                            Logger & logger)
{
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
{
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
// C - Tim, Ambient serial access = race conditions !

//    Original, 
//    for(typename Index<SymmGroup>::basis_iterator it = i.basis_begin(); !it.end(); ++it)
//        ret(0,*it,*it) = 1;

   for(std::size_t k(0); k < ret[0].n_blocks(); ++k)
       detail::iterable_matrix_impl<Matrix, SymmGroup>::left_right_boundary_init_impl(ret[0][k]);       

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
        detail::iterable_matrix_impl<Matrix, SymmGroup>::left_right_boundary_init_impl(ret[0][k]);       
    
    return ret;
}

#ifdef HAVE_ALPS_HDF5

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::load(alps::hdf5::archive & ar)
{
    ar >> alps::make_pvp("MPS",
                         static_cast<std::vector<MPSTensor<Matrix, SymmGroup> >&>(*this));
}

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::save(alps::hdf5::archive & ar) const
{
    ar << alps::make_pvp("MPS",
                         static_cast<std::vector<MPSTensor<Matrix, SymmGroup> > const &>(*this));
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
            std::cerr << "Problem on site " << i << ":" << e.what() << std::endl;
            exit(1);
        }
}


