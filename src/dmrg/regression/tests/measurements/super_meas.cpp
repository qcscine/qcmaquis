#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>
#include <iterator>
#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/block_matrix/detail/alps.hpp"

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/state_mps.h"
#include "dmrg/mp_tensors/mps_initializers.h"
#include "dmrg/mp_tensors/coherent_init.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/super_mpo.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/models/generate_mpo.hpp"

#include <boost/math/special_functions/factorials.hpp>

typedef alps::numeric::matrix<double> matrix;
typedef TrivialGroup SymmGroup;
typedef SymmGroup::charge charge;
typedef boost::tuple<charge, size_t> local_state;


std::vector<double> measure_local(MPS<matrix, SymmGroup> const& mps,
                   block_matrix<matrix, SymmGroup> const& ident, block_matrix<matrix, SymmGroup> const& op)
{
    std::vector<double> vals(mps.size());
    for (int p=0; p<mps.size(); ++p) {
        generate_mpo::MPOMaker<matrix, SymmGroup> mpom(mps.length(), ident);
        generate_mpo::Operator_Term<matrix, SymmGroup> term;
        term.operators.push_back( std::make_pair(p, op) );
        term.fill_operator = ident;
        mpom.add_term(term);
        MPO<matrix, SymmGroup> mpo = mpom.create_mpo();
        MPS<matrix, SymmGroup> super_mpo = mpo_to_smps(mpo, ident.left_basis());
        
        vals[p] = overlap(super_mpo, mps);
    }
    return vals;
}


BOOST_AUTO_TEST_CASE( density_trivial_init )
{
    maquis::cout << std::endl << "** TESTING density_trivial_init **" << std::endl;
    int L = 7;
    
    /// Bosons with Nmax=2
    const int Nmax = 2;
    SymmGroup::charge C = SymmGroup::IdentityCharge;
    Index<SymmGroup> phys_psi;
    phys_psi.insert(std::make_pair(C, Nmax+1));
    Index<SymmGroup> phys_rho = phys_psi * adjoin(phys_psi);
    
    /*
     * phys_rho = [
     *              0 --> 0,
     *              0 --> 1,
     *              0 --> 2,
     *              1 --> 0,
     *              1 --> 1,
     *              1 --> 2,
     *              2 --> 0,
     *              2 --> 1,
     *              2 --> 2,
     *            ]
     */
    
    /// building state
    std::vector<local_state> state(L);
    state[0] = local_state(C, 0);
    state[1] = local_state(C, 4);
    state[2] = local_state(C, 4);
    state[3] = local_state(C, 8);
    state[4] = local_state(C, 4);
    state[5] = local_state(C, 4);
    state[6] = local_state(C, 0);

    /// corresponding density
    std::vector<double> dens(L, 0.);
    dens[1] = 1.;
    dens[2] = 1.;
    dens[3] = 2.;
    dens[4] = 1.;
    dens[5] = 1.;

    MPS<matrix,SymmGroup> mps = state_mps<matrix>(state, phys_rho);
    
    /// operators for meas
    block_matrix<matrix, SymmGroup> ident = identity_matrix<matrix>(phys_psi);
    block_matrix<matrix, SymmGroup> densop;
    {
        matrix tmp(Nmax+1, Nmax+1, 0.);
        for (int i=1; i<Nmax+1; ++i) tmp(i,i) = i;
        densop.insert_block(tmp, C,C);
    }
    
    /// meas
    std::vector<double> meas_dens = measure_local(mps, ident, densop);
    for (int p=0; p<L; ++p) {
        maquis::cout << "site " << p << ": " << meas_dens[p] << std::endl;
        BOOST_CHECK_CLOSE(dens[p], meas_dens[p], 1e-8 );
    }
}

BOOST_AUTO_TEST_CASE( density_join_init )
{
    maquis::cout << std::endl << "** TESTING density_join_init **" << std::endl;
    int L = 7;
    
    // Bosons with Nmax=2
    const int Nmax = 2;
    SymmGroup::charge C = SymmGroup::IdentityCharge;
    Index<SymmGroup> phys_psi;
    phys_psi.insert(std::make_pair(C, Nmax+1));
    Index<SymmGroup> phys_rho = phys_psi * adjoin(phys_psi);
    
    MPS<matrix,SymmGroup> mps;
    {
    std::vector<local_state> state(L);
    state[0] = local_state(C, 4);
    state[1] = local_state(C, 4);
    state[2] = local_state(C, 4);
    state[3] = local_state(C, 4);
    state[4] = local_state(C, 4);
    state[5] = local_state(C, 4);
    state[6] = local_state(C, 4);
        
        mps =state_mps<matrix>(state, phys_rho);
    }
    
    {
        std::vector<local_state> state(L);
        state[0] = boost::make_tuple(C, 0);
        state[1] = boost::make_tuple(C, 0);
        state[2] = boost::make_tuple(C, 4);
        state[3] = boost::make_tuple(C, 4);
        state[4] = boost::make_tuple(C, 4);
        state[5] = boost::make_tuple(C, 0);
        state[6] = boost::make_tuple(C, 0);
        
        MPS<matrix,SymmGroup> tmp =state_mps<matrix>(state, phys_rho);
        mps = join(tmp, mps);
    }

    /// operators for meas
    block_matrix<matrix, SymmGroup> ident = identity_matrix<matrix>(phys_psi);
    block_matrix<matrix, SymmGroup> densop;
    {
        matrix tmp(Nmax+1, Nmax+1, 0.);
        for (int i=1; i<Nmax+1; ++i) tmp(i,i) = i;
        densop.insert_block(tmp, C,C);
    }
    
    /// meas
    std::vector<double> meas_dens = measure_local(mps, ident, densop);
    for (int p=0; p<L; ++p) {
        maquis::cout << "site " << p << ": " << meas_dens[p] << std::endl;
    }
}

BOOST_AUTO_TEST_CASE( density_coherent_init )
{
    maquis::cout << std::endl << "** TESTING density_coherent_init **" << std::endl;
    
    using std::sqrt;
    
    int L = 2;
    const int Nmax = 2;
    
    SymmGroup::charge C = SymmGroup::IdentityCharge;
    Index<SymmGroup> phys_psi;
    phys_psi.insert(std::make_pair(C, Nmax+1));
    Index<SymmGroup> phys_rho = phys_psi * adjoin(phys_psi);
    
    /// desired density
    std::vector<double> coeff(L, 0.);
    coeff[0] = sqrt(0.0075);
    coeff[1] = sqrt(0.0025);
    
    /// build coherent init MPS
    MPS<matrix,SymmGroup> mps = coherent_init_dm<matrix>(coeff, phys_psi, phys_rho);
    
    double nn = sqrt(norm(mps));
    std::cout << "norm = " << nn << std::endl;
    
    /// operators for meas
    block_matrix<matrix, SymmGroup> ident = identity_matrix<matrix>(phys_psi);
    block_matrix<matrix, SymmGroup> densop;
    {
        matrix tmp(Nmax+1, Nmax+1, 0.);
        for (int i=1; i<Nmax+1; ++i) tmp(i,i) = i;
        densop.insert_block(tmp, C,C);
    }
    
    /// meas
    std::vector<double> meas_dens = measure_local(mps, ident, densop);
    for (int p=0; p<L; ++p) {
        maquis::cout << "site " << p << ": " << meas_dens[p]/nn << std::endl;
        BOOST_CHECK_CLOSE(coeff[p]*coeff[p], meas_dens[p]/nn, 5. );
    }
}


BOOST_AUTO_TEST_CASE( density_random_init )
{
    maquis::cout << std::endl << "** TESTING density_random_init **" << std::endl;
    
    int L = 6;
    int M = 20;
    
    // Bosons with Nmax=2
    const int Nmax = 2;
    SymmGroup::charge C = SymmGroup::IdentityCharge;
    Index<SymmGroup> phys_psi;
    phys_psi.insert(std::make_pair(C, 3));
    Index<SymmGroup> phys_rho = phys_psi * adjoin(phys_psi);
        
    /// building random state
    default_mps_init<matrix, SymmGroup> initializer;
    
    MPS<matrix,SymmGroup> mps;
    mps.resize(L); initializer(mps, M, phys_rho, C);
    mps.normalize_left();
    
    /// operators for meas
    block_matrix<matrix, SymmGroup> ident = identity_matrix<matrix>(phys_psi);
    block_matrix<matrix, SymmGroup> densop;
    {
        matrix tmp(Nmax+1, Nmax+1, 0.);
        for (int i=1; i<Nmax+1; ++i) tmp(i,i) = i;
        densop.insert_block(tmp, C,C);
    }
    
    /// meas
    std::vector<double> meas_dens = measure_local(mps, densop, ident);
    for (int p=0; p<L; ++p)
        maquis::cout << "site " << p << ": " << meas_dens[p] << std::endl;
}



