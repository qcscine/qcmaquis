#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <iterator>
#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/block_matrix/detail/alps.hpp"

#include "dmrg/utils/DmrgParameters2.h"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_initializers.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/models/generate_mpo.hpp"


typedef U1 SymmGroup;
typedef alps::numeric::matrix<double> matrix;

std::ostream& operator<< (std::ostream& os, std::vector<double> const& v){
    os << "[";
    std::copy(v.begin(),v.end(),std::ostream_iterator<double>(os,"  "));
    os << "]";
    return os;
}

std::vector<double> measure_local(MPS<matrix, SymmGroup> & mps, block_matrix<matrix, SymmGroup> const & op)
{
    std::vector<double> meas(mps.length());
    for (size_t p=0; p<mps.length(); ++p) {
        mps.canonize(p);
        meas[p] = mps[p].scalar_overlap(contraction::local_op(mps[p], op));
    }
    return meas;
}

template<class Matrix, class SymmGroup>
Boundary<Matrix, SymmGroup>
mixed_left_boundary(MPS<Matrix, SymmGroup> const & bra, MPS<Matrix, SymmGroup> const & ket)
{
    assert(ket.length() == bra.length());
    Index<SymmGroup> i = ket[0].row_dim();
    Index<SymmGroup> j = bra[0].row_dim();
    Boundary<Matrix, SymmGroup> ret(i, j, 1);
    
    for(typename Index<SymmGroup>::basis_iterator it1 = i.basis_begin(); !it1.end(); ++it1)
    	for(typename Index<SymmGroup>::basis_iterator it2 = j.basis_begin(); !it2.end(); ++it2)
            ret(0, *it1, *it2) = 1;
    
    return ret;
}


template<class Matrix, class SymmGroup>
typename Matrix::value_type expval(MPS<Matrix, SymmGroup> const & bra,
                                   MPS<Matrix, SymmGroup> const & ket,
                                   MPO<Matrix, SymmGroup> const & mpo,
                                   bool verbose = false)
{
    assert(mpo.length() == bra.length() && bra.length() == ket.length());
    std::size_t L = bra.length();
    
    Boundary<Matrix, SymmGroup> left = mixed_left_boundary(bra, ket);
    
    for (int i = 0; i < L; ++i) {
        if (verbose)
            std::cout << "expval site " << i << std::endl;
        left = contraction::overlap_mpo_left_step(bra[i], ket[i], left, mpo[i]);
    }
    
    return left.traces()[0];
}

std::vector<double> measure_local(MPS<matrix, SymmGroup> const& bra, MPS<matrix, SymmGroup> const& ket,
                                  block_matrix<matrix, SymmGroup> const & op, block_matrix<matrix, SymmGroup> const & ident)
{
    assert(bra.length() == ket.length());
    
    std::vector<double> meas(bra.length());
    for (size_t p=0; p<bra.length(); ++p) {
        generate_mpo::MPOMaker<matrix, SymmGroup> mpom(bra.length(), ident);
        generate_mpo::Operator_Term<matrix, SymmGroup> term;
        term.operators.push_back( std::make_pair(p, op) );
        term.fill_operator = ident;
        mpom.add_term(term);
        MPO<matrix, SymmGroup> mpo = mpom.create_mpo();
        
        meas[p] = expval(bra, ket, mpo);
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
                     block_matrix<matrix, SymmGroup> const & ident,
                     block_matrix<matrix, SymmGroup> const & densop,
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
    assert(b1.size() == b2.size());
    int L = b1.size();
    
    // Bosons with Nmax=2
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(2, 1));
    SymmGroup::charge initc = L/2;
    
    block_matrix<matrix, SymmGroup> densop;
    densop.insert_block(matrix(1,1,1), 1,1);
    densop.insert_block(matrix(1,1,2), 2,2);
    
    block_matrix<matrix, SymmGroup> ident = identity_matrix<matrix>(phys);
    
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
    std::cout << std::endl << std::endl << "*** join_semirnd_mps_cmp_dens ***" << std::endl;
    int M = 10;
    int L = 4;
    
    std::vector<std::pair<SymmGroup::charge, size_t> > b(L, std::pair<SymmGroup::charge, size_t>(0,0));
    b[0].first = 1;
    b[2].first = 1;
    
    // Bosons with Nmax=2
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(2, 1));
    SymmGroup::charge initc = L/2;
    
    block_matrix<matrix, SymmGroup> densop;
    densop.insert_block(matrix(1,1,1), 1,1);
    densop.insert_block(matrix(1,1,2), 2,2);
    
    block_matrix<matrix, SymmGroup> ident = identity_matrix<matrix>(phys);
    
    default_mps_init<matrix, SymmGroup> initializer;

    MPS<matrix,SymmGroup> mps1;
    mps1.resize(L); initializer(mps1, M, phys, initc);
    MPS<matrix,SymmGroup> mps2 = state_mps(b, phys);
    
    run_test_bosons(L, phys, ident, densop, mps1, mps2, false);
}

BOOST_AUTO_TEST_CASE( join_rnd_mps_cmp_dens )
{
    std::cout << std::endl << std::endl << "*** join_rnd_mps_cmp_dens ***" << std::endl;
    
    int Nmps = 4;
    int M = 10;
    int L = 4;
    
    // Bosons with Nmax=2
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(2, 1));
    SymmGroup::charge initc = L/2;
    
    block_matrix<matrix, SymmGroup> densop;
    densop.insert_block(matrix(1,1,1), 1,1);
    densop.insert_block(matrix(1,1,2), 2,2);
    
    block_matrix<matrix, SymmGroup> ident = identity_matrix<matrix>(phys);
    
    default_mps_init<matrix, SymmGroup> initializer;
    
    MPS<matrix,SymmGroup> mps;
    mps.resize(L); initializer(mps, M, phys, initc);
    std::cout << "norm(mps): " << norm(mps) << std::endl;
    mps.normalize_left();
    double n = norm(mps);
    std::cout << "norm(mps) - normalized: " << n << std::endl;
    
    std::vector<double> dens = measure_local(mps, densop);
    std::cout << "Dens: " << dens << std::endl;
    
    
    for (int i=0; i<Nmps; ++i) {
        
        MPS<matrix,SymmGroup> mps2;
        mps2.resize(L); initializer(mps2, M, phys, initc);
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
