/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>

#include <iterator>
#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

#include "utils/fpcomparison.h"

#include "dmrg/block_matrix/detail/alps.hpp"

#include "dmrg/utils/DmrgParameters.h"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_initializers.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/lattice/ChainLattice.hpp"

typedef TwoU1PG SymmGroup;
typedef alps::numeric::matrix<double> matrix;

std::ostream& operator<< (std::ostream& os, std::vector<double> const& v){
    os << "[";
    std::copy(v.begin(),v.end(),std::ostream_iterator<double>(os,"  "));
    os << "]";
    return os;
}

std::vector<double> measure_local(MPS<matrix, SymmGroup> & mps, operator_selector<matrix, SymmGroup>::type const & op)
{
    std::vector<double> meas(mps.length());
    for (size_t p=0; p<mps.length(); ++p) {
        mps.canonize(p);
        meas[p] = mps[p].scalar_overlap(contraction::local_op(mps[p], op));
    }
    return meas;
}

std::vector<double> measure_local(MPS<matrix, SymmGroup> const& bra, MPS<matrix, SymmGroup> const& ket,
                                  operator_selector<matrix, SymmGroup>::type const & op,
                                  operator_selector<matrix, SymmGroup>::type const & ident)
{
    assert(bra.length() == ket.length());
    typedef operator_selector<matrix, SymmGroup>::type op_t;
    typedef std::vector<op_t> op_vec;
    std::shared_ptr<lattice_impl> lat_ptr(new ChainLattice(bra.length()));
    Lattice lattice(lat_ptr);

    std::vector<double> meas(bra.length());
    for (size_t p=0; p<bra.length(); ++p) {
        generate_mpo::MPOMaker<matrix, SymmGroup> mpom(lattice, op_vec(1,ident),op_vec(1,ident));
        generate_mpo::OperatorTerm<matrix, SymmGroup> term;
        term.operators.push_back( std::make_pair(p, op) );
        term.fill_operator = ident;
        mpom.add_term(term);
        MPO<matrix, SymmGroup> mpo = mpom.create_mpo();

        meas[p] = maquis::real(expval(bra, ket, mpo));
    }
    return meas;
}


MPS<matrix, SymmGroup> state_mps(std::vector<std::pair<SymmGroup::charge, size_t> > const & state, Index<SymmGroup> const & phys)
{
    MPS<matrix, SymmGroup> mps(state.size());

    Index<SymmGroup> curr_i;
    curr_i.insert(std::make_pair(0, 1));
    size_t curr = 0;
    for (int i=0; i<state.size(); ++i)
    {
        SymmGroup::charge newc = SymmGroup::fuse(curr_i[0].first, state[i].first);
        size_t news = (i == state.size()-1) ? 1 : phys[phys.position(state[i].first)].second;
        Index<SymmGroup> new_i;
        new_i.insert(std::make_pair(newc, news));
        ProductBasis<SymmGroup> left(phys, curr_i);
        mps[i] = MPSTensor<matrix, SymmGroup>(phys, curr_i, new_i, false, 0);
        size_t in = left(state[i].first, curr_i[0].first) + state[i].second * curr_i[0].second + curr;
        size_t out = (i == state.size()-1) ? 0 : state[i].second;
        mps[i].make_left_paired();
        mps[i].data()(SymmGroup::fuse(curr_i[0].first, state[i].first), new_i[0].first)(in, out) = 1.;
        curr_i = new_i;
        curr = state[i].second;
    }
    mps.normalize_left();
    return mps;
}

void run_test_bosons(int L, Index<SymmGroup> const& phys,
                     operator_selector<matrix, SymmGroup>::type const & ident,
                     operator_selector<matrix, SymmGroup>::type const & densop,
                     MPS<matrix,SymmGroup> & mps1, MPS<matrix,SymmGroup> & mps2,
                     bool verbose=false)
{
    std::cout << "norm(mps1): " << norm(mps1) << std::endl;
    mps1.normalize_left();
    double n1 = norm(mps1);
    std::cout << "norm(mps1) - normalized: " << n1 << std::endl;

    if (verbose) {
        std::cout << "** MPS1 **" << std::endl;
        for (int p=0; p<L; ++p)
        {
            std::cout << "SITE " << p << ":" << std::endl;
            std::cout << mps1[p] << std::endl;
        }
    }

    std::vector<double> dens1 = measure_local(mps1, densop);
    std::cout << "Dens1: " << dens1 << std::endl;


    std::cout << "norm(mps2): " << norm(mps2) << std::endl;
    mps2.normalize_left();
    double n2 = norm(mps2);
    std::cout << "norm(mps2) - normalized: " << n2 << std::endl;

    std::vector<double> dens2 = measure_local(mps2, densop);
    std::cout << "Dens2: " << dens2 << std::endl;


    double n12 = overlap(mps2, mps1);
    std::cout << "Overlap: " << n12 << std::endl;
    std::vector<double> densmix = measure_local(mps2, mps1, densop, ident);
    std::cout << "Dens-mix: " << densmix << std::endl;

    MPS<matrix,SymmGroup> mps3 = join(mps1, mps2);
    if (verbose) {
        std::cout << "** MPS3 **" << std::endl;
        for (int p=0; p<L; ++p)
        {
            std::cout << "SITE " << p << ":" << std::endl;
            std::cout << mps3[p] << std::endl;
        }
    }

    for (int p=0; p<L; ++p)
    {
        std::cout << "Making site "<<p<<" right_paired." << std::endl;
        block_matrix<matrix, SymmGroup> tmp;
        reshape_left_to_right<matrix>(phys, mps3[p].row_dim(), mps3[p].col_dim(), mps3[p].data(), tmp);
        assert( weak_equal(mps3[p].row_dim(),tmp.left_basis()) );
        assert(mps3[p].reasonable());
        std::cout<< "site " << p << " is reasonable." <<std::endl;
        if (verbose)
            std::cout << tmp << std::endl;
    }

    std::vector<double> dens3_n = measure_local(mps3, mps3, densop, ident);
    std::cout << "Dens3 no norm: " << dens3_n << std::endl;

    std::cout << "norm(mps3): " << norm(mps3) << std::endl;
    mps3.normalize_left();

    std::cout << "norm(mps3) - normalized: " << norm(mps3) << std::endl;

    std::vector<double> dens3 = measure_local(mps3, densop);
    std::cout << "Dens3: " << dens3 << std::endl;

    for (int p=0; p<L; ++p)
        BOOST_CHECK_CLOSE(dens3[p], (dens1[p]*n1+dens2[p]*n2+2*densmix[p])/(n1+n2+2*n12), 1e-8 );
}

void test_bosons(std::vector<std::pair<SymmGroup::charge, size_t> > const& b1,
                 std::vector<std::pair<SymmGroup::charge, size_t> > const& b2,
                 int M,
                 bool verbose = false)
{
    typedef operator_selector<matrix, SymmGroup>::type op_t;
    assert(b1.size() == b2.size());
    int L = b1.size();

    // Bosons with Nmax=2
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(2, 1));
    SymmGroup::charge initc = L/2;

    op_t densop;
    densop.insert_block(matrix(1,1,1), 1,1);
    densop.insert_block(matrix(1,1,2), 2,2);

    op_t ident = identity_matrix<op_t>(phys);

    MPS<matrix,SymmGroup> mps1 = state_mps(b1, phys);
    MPS<matrix,SymmGroup> mps2 = state_mps(b2, phys);
    run_test_bosons(L, phys, ident, densop, mps1, mps2, verbose);
}


BOOST_AUTO_TEST_CASE( join_orthogonal_mps_cmp_dens )
{
    std::cout << std::endl << std::endl << "*** join_orthogonal_mps_cmp_dens ***" << std::endl;
    int M = 10;
    int L = 4;

    std::vector<std::pair<SymmGroup::charge, size_t> > b1(L, std::pair<SymmGroup::charge, size_t>(0,0));
    b1[0].first = 1;
    b1[2].first = 1;

    std::vector<std::pair<SymmGroup::charge, size_t> > b2(L, std::pair<SymmGroup::charge, size_t>(0,0));
    b2[1].first = 1;
    b2[3].first = 1;

    test_bosons(b1, b2, M, false);
}

BOOST_AUTO_TEST_CASE( join_mps_cmp_dens )
{
    std::cout << std::endl << std::endl << "*** join_mps_cmp_dens ***" << std::endl;
    int M = 10;
    int L = 4;

    std::vector<std::pair<SymmGroup::charge, size_t> > b1(L, std::pair<SymmGroup::charge, size_t>(0,0));
    b1[0].first = 1;
    b1[2].first = 1;

    std::vector<std::pair<SymmGroup::charge, size_t> > b2(L, std::pair<SymmGroup::charge, size_t>(0,0));
    b2[0].first = 1;
    b2[3].first = 1;

    test_bosons(b1, b2, M, false);
}

BOOST_AUTO_TEST_CASE( join_same_mps_cmp_dens )
{
    std::cout << std::endl << std::endl << "*** join_same_mps_cmp_dens ***" << std::endl;
    int M = 10;
    int L = 4;

    std::vector<std::pair<SymmGroup::charge, size_t> > b1(L, std::pair<SymmGroup::charge, size_t>(0,0));
    b1[0].first = 1;
    b1[2].first = 1;

    test_bosons(b1, b1, M, false);
}


BOOST_AUTO_TEST_CASE( join_semirnd_mps_cmp_dens )
{
    typedef operator_selector<matrix, SymmGroup>::type op_t;
    std::cout << std::endl << std::endl << "*** join_semirnd_mps_cmp_dens ***" << std::endl;
    int M = 10;
    int L = 4;

    DmrgParameters parms;

    std::vector<std::pair<SymmGroup::charge, size_t> > b(L, std::pair<SymmGroup::charge, size_t>(0,0));
    b[0].first = 1;
    b[2].first = 1;

    // Bosons with Nmax=2
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(2, 1));
    SymmGroup::charge initc = L/2;

    op_t densop;
    densop.insert_block(matrix(1,1,1), 1,1);
    densop.insert_block(matrix(1,1,2), 2,2);

    op_t ident = identity_matrix<op_t>(phys);

    parms.set("max_bond_dimension", M);
    default_mps_init<matrix, SymmGroup> initializer(parms, std::vector<Index<SymmGroup> >(1, phys), initc, std::vector<int>(L,0));

    MPS<matrix,SymmGroup> mps1;
    mps1.resize(L); initializer(mps1);
    MPS<matrix,SymmGroup> mps2 = state_mps(b, phys);

    run_test_bosons(L, phys, ident, densop, mps1, mps2, false);
}

BOOST_AUTO_TEST_CASE( join_rnd_mps_cmp_dens )
{
    typedef operator_selector<matrix, SymmGroup>::type op_t;
    std::cout << std::endl << std::endl << "*** join_rnd_mps_cmp_dens ***" << std::endl;

    int Nmps = 4;
    int M = 10;
    int L = 4;

    DmrgParameters parms;

    // Bosons with Nmax=2
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(2, 1));
    SymmGroup::charge initc = L/2;

    op_t densop;
    densop.insert_block(matrix(1,1,1), 1,1);
    densop.insert_block(matrix(1,1,2), 2,2);

    op_t ident = identity_matrix<op_t>(phys);

    parms.set("max_bond_dimension", M);
    default_mps_init<matrix, SymmGroup> initializer(parms, std::vector<Index<SymmGroup> >(1, phys), initc, std::vector<int>(L,0));

    MPS<matrix,SymmGroup> mps;
    mps.resize(L); initializer(mps);
    std::cout << "norm(mps): " << norm(mps) << std::endl;
    mps.normalize_left();
    double n = norm(mps);
    std::cout << "norm(mps) - normalized: " << n << std::endl;

    std::vector<double> dens = measure_local(mps, densop);
    std::cout << "Dens: " << dens << std::endl;


    for (int i=0; i<Nmps; ++i) {

        MPS<matrix,SymmGroup> mps2;
        mps2.resize(L); initializer(mps2);
        std::cout << "norm(mps2): " << norm(mps2) << std::endl;
        mps2.normalize_left();
        double n2 = norm(mps2);
        std::cout << "norm(mps2) - normalized: " << n2 << std::endl;
        std::vector<double> dens2 = measure_local(mps2, densop);
        std::cout << "Dens2: " << dens2 << std::endl;

        std::vector<double> densmix = measure_local(mps2, mps, densop, ident);
        double n12 = overlap(mps,mps2);
        std::cout << "Overlap=" << n12 << std::endl;
        std::cout << "Dens-mix: " << densmix << std::endl;

        MPS<matrix,SymmGroup> mps3 = join(mps, mps2);
        std::cout << mps3.description() << std::endl;
        for (int p=0; p<mps3.length(); ++p)
            assert(mps3[p].reasonable());
        std::cout << "norm(mps3): " << norm(mps3) << std::endl;
        mps3.normalize_left();
        for (int p=0; p<mps3.length(); ++p)
            assert(mps3[p].reasonable());
        double n3 = norm(mps3);
        std::cout << "norm(mps3) - normalized: " << n3 << std::endl;
        std::vector<double> dens3 = measure_local(mps3, densop);
        std::cout << "Dens3: " << dens3 << std::endl;

        for (int p=0; p<L; ++p)
            BOOST_CHECK_CLOSE(dens3[p], (dens[p]*n+dens2[p]*n2+2*densmix[p])/(n+n2+2*n12), 1e-8 );
        std::swap(mps,mps3);
        std::swap(dens,dens3);
        std::swap(n,n3);
    }

}
