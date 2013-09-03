/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

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
    //     it turned out, this is also quite dangerous: if a block is 1x2,
    //     normalize_left will resize it to 1x1
    //     init() should take care of it, in case needed. Otherwise some
    //     adhoc states will be broken (e.g. identity MPS)
    // for (int i = 0; i < L; ++i)
    //     (*this)[i].normalize_left(DefaultSolver());

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
    #ifdef AMBIENT
    locale::compact(length()); locale l(length()-1,length()-1);
    #endif
    block_matrix<Matrix, SymmGroup> t = (*this)[length()-1].normalize_left(DefaultSolver());
    // now state is: A A A A A A A
    canonized_i = length()-1;
}

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::normalize_right()
{
    canonize(0);
    // now state is: M B B B B B B
    #ifdef AMBIENT
    locale::compact(length()); locale l(0,0);
    #endif
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
        #ifdef AMBIENT
        locale::compact(length()); locale l(i,i);
        #endif
        block_matrix<Matrix, SymmGroup> t = (*this)[i].normalize_left(method);
        if (i < length()-1) { 
             #ifdef AMBIENT
             ++l;
             #endif
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
        #ifdef AMBIENT
        locale::compact(length()); locale l(i,i);
        #endif
        block_matrix<Matrix, SymmGroup> t = (*this)[i].normalize_right(method);
        if (i > 0) { 
            #ifdef AMBIENT
            --l;
            #endif
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
template<class OtherMatrix>
truncation_results
MPS<Matrix, SymmGroup>::grow_l2r_sweep(MPOTensor<Matrix, SymmGroup> const & mpo,
                                       Boundary<OtherMatrix, SymmGroup> const & left,
                                       Boundary<OtherMatrix, SymmGroup> const & right,
                                       std::size_t l, double alpha,
                                       double cutoff, std::size_t Mmax)
{ // canonized_i invalided through (*this)[]
    MPSTensor<Matrix, SymmGroup> new_mps;
    truncation_results trunc;
    
    boost::tie(new_mps, trunc) =
    contraction::predict_new_state_l2r_sweep((*this)[l], mpo, left, right, alpha, cutoff, Mmax);
    
    (*this)[l+1] = contraction::predict_lanczos_l2r_sweep((*this)[l+1],
                                                          (*this)[l], new_mps);
    (*this)[l] = new_mps;
    return trunc;
}

template<class Matrix, class SymmGroup>
template<class OtherMatrix>
truncation_results
MPS<Matrix, SymmGroup>::grow_r2l_sweep(MPOTensor<Matrix, SymmGroup> const & mpo,
                                       Boundary<OtherMatrix, SymmGroup> const & left,
                                       Boundary<OtherMatrix, SymmGroup> const & right,
                                       std::size_t l, double alpha,
                                       double cutoff, std::size_t Mmax)
{ // canonized_i invalided through (*this)[]
    MPSTensor<Matrix, SymmGroup> new_mps;
    truncation_results trunc;
    
    boost::tie(new_mps, trunc) =
    contraction::predict_new_state_r2l_sweep((*this)[l], mpo, left, right, alpha, cutoff, Mmax);
    
    (*this)[l-1] = contraction::predict_lanczos_r2l_sweep((*this)[l-1],
                                                          (*this)[l], new_mps);
    (*this)[l] = new_mps;
    return trunc;
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

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::apply(block_matrix<Matrix, SymmGroup> const& op, MPS<Matrix, SymmGroup>::size_type p)
{
    typedef typename SymmGroup::charge charge;
    using std::size_t;
    
    /// Compute (and check) charge difference
    charge diff = SymmGroup::IdentityCharge;
    if (op.n_blocks() > 0)
        diff = SymmGroup::fuse(op.right_basis()[0].first, -op.left_basis()[0].first);
    for (size_t n=0; n< op.n_blocks(); ++n) {
        if ( SymmGroup::fuse(op.right_basis()[n].first, -op.left_basis()[n].first) != diff )
            throw std::runtime_error("Operator not allowed. All non-zero blocks have to provide same `diff`.");
    }
    
    /// Apply operator
    (*this)[p] = contraction::multiply_with_op((*this)[p], op);
    
    /// Propagate charge difference
    for (size_t i=p+1; i<length(); ++i) {
        (*this)[i].shift_aux_charges(diff);
    }
}

template<class Matrix, class SymmGroup>
void MPS<Matrix, SymmGroup>::apply(block_matrix<Matrix, SymmGroup> const& fill, block_matrix<Matrix, SymmGroup> const& op, MPS<Matrix, SymmGroup>::size_type p)
{
    for (size_t i=0; i<p; ++i) {
        (*this)[i] = contraction::multiply_with_op((*this)[i], fill);
    }
    apply(op, p);
}

template<class Matrix, class SymmGroup>
void load(std::string const& dirname, MPS<Matrix, SymmGroup> & mps)
{
    /// get size of MPS
    std::size_t L = 0;
    while (boost::filesystem::exists( dirname + "/mps" + boost::lexical_cast<std::string>(++L) + ".h5" ));
    
    /// load tensors
    MPS<Matrix, SymmGroup> tmp(L);
    size_t loop_max = tmp.length();
    semi_parallel_for(locale::compact(loop_max), locale k = 0; k < loop_max; ++k){
        std::string fname = dirname+"/mps"+boost::lexical_cast<std::string>((size_t)k)+".h5";
        storage::archive ar(fname);
        ar["/tensor"] >> tmp[k];
    }
    swap(mps, tmp);
}

template<class Matrix, class SymmGroup>
void save(std::string const& dirname, MPS<Matrix, SymmGroup> const& mps)
{
    size_t loop_max = mps.length();
#ifdef AMBIENT
    semi_parallel_for(locale::compact(loop_max), locale k = 0; k < loop_max; ++k)
        mps[k].make_left_paired();
    ambient::sync();
#endif
    semi_parallel_for(locale::compact(loop_max), locale k = 0; k < loop_max; ++k){
#ifdef AMBIENT
        if(!ambient::controller.local()) continue;
#endif
        std::string fname = dirname+"/mps"+boost::lexical_cast<std::string>((size_t)k)+".h5";
        storage::archive ar(fname, "w");
        ar["/tensor"] << mps[k];
    }
}

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


