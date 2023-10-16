/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MPS_MPO_OPS_H
#define MPS_MPO_OPS_H

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/twositetensor.h"
#include "dmrg/mp_tensors/special_mpos.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/utils/utils.hpp"
#include "utils/traits.hpp"
#include "mps_mpo_detail.h"

// Forward declaration
template<class Matrix, class SymmGroup>
typename Matrix::value_type expval(MPS<Matrix, SymmGroup> const & bra, MPS<Matrix, SymmGroup> const & ket,
                                   MPO<Matrix, SymmGroup> const & mpo);

template<class Matrix, class SymmGroup>
typename Matrix::value_type expvalFromRight(MPS<Matrix, SymmGroup> const & bra, MPS<Matrix, SymmGroup> const & ket,
                                            MPO<Matrix, SymmGroup> const & mpo);

/**
 * @brief Method to calculate the expectation value choosing the starting site.
 *
 * The expectation value is calculated by contracting the full MPS/MPO network.
 * This can be done starting either from the first or from the last site of the
 * lattice. The first option is activated with d == 0, the second one with
 * d != 0
 *
 * @param mps Input MPS
 * @param mpo Input MPO
 * @param d If == 0, starts from the first site, otherwise, starts from the last site.
 * @return Matrix::value_type Expectation value < mps | H | mps >
 */
template<class Matrix, class SymmGroup>
auto expval(MPS<Matrix, SymmGroup> const & mps, MPO<Matrix, SymmGroup> const & mpo, int d)
{
    return (d == 0) ? expval(mps, mps, mpo) : expvalFromRight(mps, mps, mpo);
}

/**
 * @brief Method to calculate the matrix element of an operator between two different MPSs
 *
 * @param bra Input MPS representing the bra
 * @param ket Input MPS representing the ket
 * @param mpo Input MPO
 * @return Matrix::value_type  < bra | mpo | ket >
 */
template<class Matrix, class SymmGroup>
typename Matrix::value_type expval(MPS<Matrix, SymmGroup> const & bra, MPS<Matrix, SymmGroup> const & ket,
                                   MPO<Matrix, SymmGroup> const & mpo)
{
    parallel::scheduler_balanced scheduler(bra.length());
    assert(mpo.length() == bra.length() && bra.length() == ket.length());
    std::size_t L = bra.length();
    Boundary<Matrix, SymmGroup> left = mps_mpo_detail::mixed_left_boundary(bra, ket);
    for (int i = 0; i < L; ++i) {
        parallel::guard proc(scheduler(i));
        left = contraction::Engine<Matrix, Matrix, SymmGroup>::overlap_mpo_left_step(bra[i], ket[i], left, mpo[i], false);
    }
    return left.traces()[0] + mpo.getCoreEnergy()*overlap(bra, ket);
}

/**
 * @brief Method to calculate the matrix element of an operator between two different MPSs
 *
 * Unlike [expVal], this method contracts the network starting from the last site.
 *
 * @param bra Input MPS representing the bra
 * @param ket Input MPS representing the ket
 * @param mpo Input MPO
 * @return Matrix::value_type  < bra | mpo | ket >
 */
template<class Matrix, class SymmGroup>
typename Matrix::value_type expvalFromRight(MPS<Matrix, SymmGroup> const & bra, MPS<Matrix, SymmGroup> const & ket,
                                            MPO<Matrix, SymmGroup> const & mpo)
{
    parallel::scheduler_balanced scheduler(bra.length());
    assert(mpo.length() == bra.length() && bra.length() == ket.length());
    std::size_t L = bra.length();
    Boundary<Matrix, SymmGroup> right = mps_mpo_detail::mixed_right_boundary(bra, ket);
    for (int i = L-1; i >= 0; --i) {
        parallel::guard proc(scheduler(i));
        right = contraction::Engine<Matrix, Matrix, SymmGroup>::overlap_mpo_right_step(bra[i], ket[i], right, mpo[i], false);
    }
    return right.traces()[0] + mpo.getCoreEnergy()*overlap(bra, ket);
}

/**
 * @brief Calculates the expectation value of an MPS over an MPO.
 * @param mps Input MPS
 * @param mpo Input MPO
 * @return Matrix::value_type < mps | mpo | mps >
 */
template<class Matrix, class SymmGroup>
typename Matrix::value_type expval(MPS<Matrix, SymmGroup> const & mps, MPO<Matrix, SymmGroup> const & mpo)
{
    return expval(mps, mps, mpo);
}

template<class Matrix, class SymmGroup>
std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> multi_expval(MPS<Matrix, SymmGroup> const & bra,
                                                                       MPS<Matrix, SymmGroup> const & ket,
                                                                       MPO<Matrix, SymmGroup> const & mpo)
{
    assert(bra.length() == ket.length());
    assert(mpo.length() == bra.length());
    std::size_t L = bra.length();

    //Boundary<Matrix, SymmGroup> left = make_left_boundary(bra, ket);
    Boundary<Matrix, SymmGroup> left = mps_mpo_detail::mixed_left_boundary(bra, ket);

    for (int i = 0; i < L; ++i)
        left = contraction::Engine<Matrix, Matrix, SymmGroup>::overlap_mpo_left_step(bra[i], ket[i], left, mpo[i]);

    return left.traces();
}

template<class Matrix, class SymmGroup>
std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> multi_expval(MPS<Matrix, SymmGroup> const & mps,
                                                                       MPO<Matrix, SymmGroup> const & mpo)
{
    return multi_expval(mps, mps, mpo);
}

template<class Matrix, class SymmGroup>
double norm(MPS<Matrix, SymmGroup> const & mps)
{
    parallel::scheduler_balanced scheduler(mps.length());
    std::size_t L = mps.length();

    block_matrix<Matrix, SymmGroup> left;
    left.insert_block(Matrix(1, 1, 1), SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);

    for(size_t i = 0; i < L; ++i) {
        parallel::guard proc(scheduler(i));
        MPSTensor<Matrix, SymmGroup> cpy = mps[i];
        left = contraction::Engine<Matrix, Matrix, SymmGroup>::overlap_left_step(mps[i], cpy, left); // serial
    }

    return maquis::real(trace(left));
}

template<class Matrix, class SymmGroup>
typename MPS<Matrix, SymmGroup>::scalar_type overlap(MPS<Matrix, SymmGroup> const & braMps,
                                                     MPS<Matrix, SymmGroup> const & ketMps)
{
    parallel::scheduler_balanced scheduler(braMps.length());
    assert(braMps.length() == ketMps.length());
    block_matrix<Matrix, SymmGroup> left;
    left.insert_block(Matrix(1, 1, 1), SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);
    for(size_t i = 0; i < braMps.length(); ++i) {
        parallel::guard proc(scheduler(i));
        left = contraction::Engine<Matrix, Matrix, SymmGroup>::overlap_left_step(braMps[i], ketMps[i], left);
    }
    return trace(left);
}

template<class Matrix, class SymmGroup>
std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> multi_overlap(MPS<Matrix, SymmGroup> const & braMps,
                                                                        MPS<Matrix, SymmGroup> const & ketMps)
{
    assert(braMps.length() == ketMps.length());
    block_matrix<Matrix, SymmGroup> left;
    left.insert_block(Matrix(1, 1, 1), SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);
    for (int i = 0; i < braMps.length(); ++i)
        left = contraction::Engine<Matrix, Matrix, SymmGroup>::overlap_left_step(braMps[i], ketMps[i], left);
    assert(left.right_basis().sum_of_sizes() == 1);
    std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> vals;
    vals.reserve(left.basis().sum_of_left_sizes());
    for (int n=0; n<left.n_blocks(); ++n)
        for (int i=0; i<left.basis().left_size(n); ++i)
            vals.push_back( left[n](i,0) );
    return vals;
}

//typedef std::vector< std::vector< std::pair<std::string, double> > > entanglement_spectrum_type;
typedef std::vector< std::pair<std::vector<std::string>, std::vector<double> > > entanglement_spectrum_type;
template<class Matrix, class SymmGroup>
std::vector<double>
calculate_bond_renyi_entropies(MPS<Matrix, SymmGroup> mps, double n,
                               std::vector<int> * measure_es_where = NULL,
                               entanglement_spectrum_type * spectra = NULL) // to be optimized later
{
    std::size_t L = mps.length();
    std::vector<double> ret;

    MPS<Matrix, SymmGroup> const& constmps = mps;

    block_matrix<Matrix, SymmGroup> lb;

    if (spectra != NULL)
        spectra->clear();

    mps.canonize(0);
    for (std::size_t p = 1; p < L; ++p)
    {
        block_matrix<Matrix, SymmGroup> t, u, v;
        block_matrix<typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type, SymmGroup> s;

        constmps[p-1].make_left_paired();
        constmps[p].make_right_paired();

        gemm(constmps[p-1].data(), constmps[p].data(), t);

        svd(t, u, v, s);

        std::vector<double> sv = maquis::dmrg::detail::bond_renyi_entropies(s);

        if (spectra != NULL && measure_es_where != NULL
            && std::find(measure_es_where->begin(), measure_es_where->end(), p) != measure_es_where->end()) {
            std::vector< std::string > labels;
            std::vector< double > values;
            for (std::size_t k = 0; k < s.n_blocks(); ++k) {
                std::ostringstream oss_c;
                oss_c << s.left_basis()[k].first;
                std::string c_str = oss_c.str();
                for (std::size_t l = 0; l < s.left_basis()[k].second; ++l) {
                    labels.push_back( c_str );
                    values.push_back( s[k](l,l) );
                }
            }
            spectra->push_back(std::make_pair(labels, values));
        }

        double S = 0;
        if (n == 1) {
            for (std::vector<double>::const_iterator it = sv.begin();
                 it != sv.end(); ++it)
                S += *it * log(*it);
            ret.push_back(-S);
        } else {
            for (std::vector<double>::const_iterator it = sv.begin();
                 it != sv.end(); ++it)
                S += pow(*it, n);
            ret.push_back(1/(1-n)*log(S));
        }

        mps.move_normalization_l2r(p-1, p, DefaultSolver());
    }

    return ret;
}

template<class Matrix, class SymmGroup>
std::vector<double>
calculate_bond_entropies(MPS<Matrix, SymmGroup> & mps)
{
    return calculate_bond_renyi_entropies(mps, 1, NULL);
}

template<class Matrix, class SymmGroup>
typename MPS<Matrix, SymmGroup>::scalar_type dm_trace(MPS<Matrix, SymmGroup> const& mps, Index<SymmGroup> const& phys_psi)
{
    typedef typename SymmGroup::charge charge;
    charge I = SymmGroup::IdentityCharge;
    size_t L = mps.length();

    Index<SymmGroup> phys_rho = phys_psi * adjoin(phys_psi);
    ProductBasis<SymmGroup> pb(phys_psi, phys_psi, boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                       boost::lambda::_1, -boost::lambda::_2));

    Matrix identblock(phys_rho.size_of_block(I), 1, 0.);
    for (int s=0; s<phys_psi.size(); ++s)
        for (int ss=0; ss<phys_psi[s].second; ++ss) {
            identblock(pb(phys_psi[s].first, phys_psi[s].first) + ss*phys_psi[s].second+ss, 0) = 1.;
        }
    block_matrix<Matrix, SymmGroup> ident;
    ident.insert_block(identblock, I, I);

    Index<SymmGroup> trivial_i;
    trivial_i.insert(std::make_pair(I, 1));
    MPSTensor<Matrix, SymmGroup> mident(phys_rho, trivial_i, trivial_i);
    mident.data() = ident;

    MPS<Matrix,SymmGroup> mps_ident(L);
    for (int p=0; p<L; ++p)
        mps_ident[p] = mident;

    return overlap(mps, mps_ident);
}


// Specific to Fermi-Hubbard on a Ladder!!
template<class Matrix, class SymmGroup>
void fix_density(MPS<Matrix, SymmGroup> & mps, std::vector<typename operator_selector<Matrix, SymmGroup>::type> const & dens_ops,
                 std::vector<std::vector<double> > const & dens)
{
    typedef typename operator_selector<Matrix, SymmGroup>::type op_t;

    assert( mps.size() == dens[0].size() );
    assert( dens_ops.size() == dens.size() );
    size_t L = mps.size();

    mps.normalize_left();
    mps.canonize(0);
    for (int p=0; p<L; ++p)
    {
        Index<SymmGroup> phys = mps[p].site_dim();
        typename SymmGroup::charge empty, up, down, updown;
        empty[0]  = 0;  empty[1]  = 0;
        up[0]     = 1;  up[1]     = 0;
        down[0]   = 0;  down[1]   = 1;
        updown[0] = 1;  updown[1] = 1;


        block_matrix<Matrix, SymmGroup> rho = contraction::density_matrix(mps[p], mps[p]);

        for (size_t j=0; j<dens.size(); ++j) {

            MPSTensor<Matrix, SymmGroup> tmp = contraction::local_op(mps[p], dens_ops[j]);
            double cur_dens = mps[p].scalar_overlap(tmp);
            maquis::cout << "Density[" << j << "] (before) = " << cur_dens << std::endl;
        }

        double a = trace(rho(down, down)) * trace(rho(updown, updown));
        double b = trace(rho(up, up)) * trace(rho(down, down)) + dens[0][p] * trace(rho(updown, updown)) - dens[1][p] * trace(rho(updown, updown));
        double c = - dens[1][p] * trace(rho(up, up));
        double k2 = ( -b + sqrt(b*b - 4*a*c) ) / (2*a);

        double k1 = dens[0][p] / ( trace(rho(up, up)) + k2*trace(rho(updown,updown)) );

        double t0 = 0.;
        t0 += k1*trace( rho(up, up) );
        t0 += k2*trace( rho(down, down) );
        t0 += k1*k2*trace( rho(updown, updown) );
        double k0 = (1.-t0) / trace(rho(empty, empty));

        maquis::cout << "k0 = " << k0 << std::endl;
        maquis::cout << "k1 = " << k1 << std::endl;
        maquis::cout << "k2 = " << k2 << std::endl;
        assert( k0 > 0 ); // not always the case!!!

        op_t rescale = identity_matrix<typename operator_selector<Matrix, SymmGroup>::type>(phys);
        rescale(empty, empty) *= std::sqrt(k0);
        rescale(up, up) *= std::sqrt(k1);
        rescale(down, down) *= std::sqrt(k2);
        rescale(updown, updown) *= std::sqrt(k1*k2);

        mps[p] = contraction::local_op(mps[p], rescale);

        {
            for (size_t j=0; j<dens.size(); ++j) {
                MPSTensor<Matrix, SymmGroup> tmp = contraction::local_op(mps[p], dens_ops[j]);
                double meas_dens = mps[p].scalar_overlap(tmp) / mps[p].scalar_norm();
                maquis::cout << "Density[" << j << "] (after) = " << meas_dens << ", should be " << dens[j][p] << std::endl;
            }
        }

        block_matrix<Matrix, SymmGroup> t_norm = mps[p].normalize_left(DefaultSolver());
        if (p < L-1)
            mps[p+1].multiply_from_left(t_norm);

    }

}


#endif
