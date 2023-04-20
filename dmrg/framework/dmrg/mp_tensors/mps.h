/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MPS_H
#define MPS_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/mp_tensors/boundary.h"

#include <limits>

template<class Matrix, class SymmGroup>
struct mps_initializer;

template<class Matrix, class SymmGroup>
class MPS
{
    typedef std::vector<MPSTensor<Matrix, SymmGroup> > data_t;
public:
    typedef std::size_t size_t;

    // reproducing interface of std::vector
    typedef typename data_t::size_type size_type;
    typedef typename data_t::value_type value_type;
    typedef typename data_t::iterator iterator;
    typedef typename data_t::const_iterator const_iterator;
    typedef typename MPSTensor<Matrix, SymmGroup>::scalar_type scalar_type;

    MPS();
    MPS(size_t L);
    MPS(size_t L, mps_initializer<Matrix, SymmGroup> & init);
    MPS(std::initializer_list<MPSTensor<Matrix, SymmGroup> > l);

    size_t size() const { return data_.size(); }
    size_t length() const { return size(); }
    Index<SymmGroup> const & site_dim(size_t i) const { return data_[i].site_dim(); }
    Index<SymmGroup> const & row_dim(size_t i) const { return data_[i].row_dim(); }
    Index<SymmGroup> const & col_dim(size_t i) const { return data_[i].col_dim(); }

    value_type const & operator[](size_t i) const;
    value_type& operator[](size_t i);

    void resize(size_t L);

    const_iterator begin() const {return data_.begin();}
    const_iterator end() const {return data_.end();}
    const_iterator const_begin() const {return data_.begin();}
    const_iterator const_end() const {return data_.end();}
    iterator begin() {return data_.begin();}
    iterator end() {return data_.end();}

    size_t canonization(bool=false) const;
    void canonize(size_t center, DecompMethod method = DefaultSolver());

    void setComplexPartToZero();
    void scaleByScalar(scalar_type scalingFactor);

    void normalize_left();
    void normalize_right();

    void move_normalization_l2r(size_t p1, size_t p2, DecompMethod method=DefaultSolver());
    void move_normalization_r2l(size_t p1, size_t p2, DecompMethod method=DefaultSolver());

    std::string description() const;

    template<class OtherMatrix>
    truncation_results grow_l2r_sweep(MPOTensor<Matrix, SymmGroup> const & mpo,
                                      Boundary<OtherMatrix, SymmGroup> const & left,
                                      Boundary<OtherMatrix, SymmGroup> const & right,
                                      std::size_t l, double alpha,
                                      double cutoff, std::size_t Mmax,
                                      bool perturbDM, bool verbose);
    template<class OtherMatrix>
    truncation_results grow_r2l_sweep(MPOTensor<Matrix, SymmGroup> const & mpo,
                                      Boundary<OtherMatrix, SymmGroup> const & left,
                                      Boundary<OtherMatrix, SymmGroup> const & right,
                                      std::size_t l, double alpha,
                                      double cutoff, std::size_t Mmax,
                                      bool perturbDM, bool verbose);

    Boundary<Matrix, SymmGroup> left_boundary() const;
    Boundary<Matrix, SymmGroup> right_boundary() const;

    void apply(typename operator_selector<Matrix, SymmGroup>::type const&, size_type);
    void apply(typename operator_selector<Matrix, SymmGroup>::type const&,
               typename operator_selector<Matrix, SymmGroup>::type const&, size_type);

    friend void swap(MPS& a, MPS& b)
    {
        using std::swap;
        swap(a.data_, b.data_);
        swap(a.canonized_i, b.canonized_i);
    }

    template <class Archive> void serialize(Archive & ar, const unsigned int version);

private:

    data_t data_;
    mutable size_t canonized_i;
};

template<class Matrix, class SymmGroup>
void load(std::string const& dirname, MPS<Matrix, SymmGroup> & mps);
template<class Matrix, class SymmGroup>
void save(std::string const& dirname, MPS<Matrix, SymmGroup> const& mps);

template<class Matrix, class SymmGroup>
struct mps_initializer
{
    virtual ~mps_initializer() {}
    virtual void operator()(MPS<Matrix, SymmGroup> & mps) = 0;
};

template<class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup> join(MPS<Matrix, SymmGroup> const & a,
                            MPS<Matrix, SymmGroup> const & b,
                            double alpha=1., double beta=1.)
{
    assert( a.length() == b.length() );

    MPSTensor<Matrix, SymmGroup> aright=a[a.length()-1], bright=b[a.length()-1];
    aright.multiply_by_scalar(alpha);
    bright.multiply_by_scalar(beta);

    MPS<Matrix, SymmGroup> ret(a.length());
    ret[0] = join(a[0],b[0],l_boundary_f);
    ret[a.length()-1] = join(aright,bright,r_boundary_f);
    for (std::size_t p = 1; p < a.length()-1; ++p)
        ret[p] = join(a[p], b[p]);
    return ret;
}

template<class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup> join_general(MPS<Matrix, SymmGroup> const & a,
                                    MPS<Matrix, SymmGroup> const & b,
                                    typename MPSTensor<Matrix, SymmGroup>::scalar_type alpha=1.,
                                    typename MPSTensor<Matrix, SymmGroup>::scalar_type beta=1.)

{
    assert( a.length() == b.length() );

    MPSTensor<Matrix, SymmGroup> aright=a[a.length()-1], bright=b[a.length()-1];
    aright.multiply_by_scalar(alpha);
    bright.multiply_by_scalar(beta);

    MPS<Matrix, SymmGroup> ret(a.length());
    ret[0] = join(a[0],b[0],l_boundary_f);
    ret[a.length()-1] = join(aright,bright,r_boundary_f);
    for (std::size_t p = 1; p < a.length()-1; ++p)
        ret[p] = join(a[p], b[p]);
    return ret;
}

#include "dmrg/mp_tensors/compression.h"

template<class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup> joinAndTruncate(MPS<Matrix, SymmGroup> & a,
                                       MPS<Matrix, SymmGroup> & b,
                                       int mMax)

{
    assert( a.length() == b.length() );
    int nOfSites = a.length();
    MPS<Matrix, SymmGroup> ret(nOfSites);
#pragma omp parallel for
    for (int p = 0; p < nOfSites; ++p) {
        if (p == 0)
            ret[0] = join(a[0], b[0], l_boundary_f);
        else if (p == nOfSites-1)
            ret[nOfSites-1] = join(a[nOfSites-1], b[nOfSites-1], r_boundary_f);
        else
            ret[p] = join(a[p], b[p]);
    }
    ret = compression::l2r_compress(ret, mMax, 0.);
    return ret;
}

template<class Matrix, class SymmGroup>
Boundary<Matrix, SymmGroup>
make_left_boundary(MPS<Matrix, SymmGroup> const & bra, MPS<Matrix, SymmGroup> const & ket)
{
    assert(ket.length() == bra.length());
    Index<SymmGroup> i = ket[0].row_dim();
    Index<SymmGroup> j = bra[0].row_dim();
    Boundary<Matrix, SymmGroup> ret(i, j, 1);

    for(typename Index<SymmGroup>::basis_iterator it1 = i.basis_begin(); !it1.end(); ++it1)
        for(typename Index<SymmGroup>::basis_iterator it2 = j.basis_begin(); !it2.end(); ++it2)
            ret[0](*it1, *it2) = 1;

    return ret;
}

template<class Matrix, class SymmGroup>
Boundary<Matrix, SymmGroup>
make_right_boundary(MPS<Matrix, SymmGroup> const & bra, MPS<Matrix, SymmGroup> const & ket)
{
    assert(ket.length() == bra.length());
    std::size_t L = ket.length();
    Index<SymmGroup> i = ket[L-1].col_dim();
    Index<SymmGroup> j = bra[L-1].col_dim();
    Boundary<Matrix, SymmGroup> ret(j, i, 1);

    for(typename Index<SymmGroup>::basis_iterator it1 = i.basis_begin(); !it1.end(); ++it1)
        for(typename Index<SymmGroup>::basis_iterator it2 = j.basis_begin(); !it2.end(); ++it2)
            ret[0](*it2, *it1) = 1;

    return ret;
}

#include "dmrg/mp_tensors/mps.hpp"

#endif
