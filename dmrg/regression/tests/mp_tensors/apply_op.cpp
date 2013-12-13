/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
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
#include <boost/test/floating_point_comparison.hpp>
#include <boost/mpl/list.hpp>

#include <iterator>
#include <iostream>


using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/block_matrix/detail/alps.hpp"

#include "dmrg/utils/DmrgParameters.h"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_initializers.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/coded/lattice.hpp"


typedef alps::numeric::matrix<double> matrix;


struct U1System {
    
    typedef U1 grp;
    
    static const bool fermionic()   { return false; }
    static const int max_bond_dim() { return 20; }
    
    static const std::string lattice_lib() { return "alps";  }
    static const std::string model_lib()   { return "coded"; }
    
    static DmrgParameters parms()
    {
        DmrgParameters p;
        p.set("max_bond_dimension", max_bond_dim());
        p.set("lattice_library", lattice_lib());
        p.set("model_library", model_lib());
        p << model_parms();
        return p;
    }
    static ModelParameters model_parms()
    {
        ModelParameters p;
        p.set("LATTICE", "open chain lattice");
        p.set("L", 10);
        
        p.set("MODEL",           "boson Hubbard");
        p.set("Nmax",            2              );
        p.set("t",               1.             );
        p.set("U",               1.             );
        p.set("u1_total_charge", 4              );
        
        return p;
    }
    
    static const std::vector<std::string> dens_ops()
    {
        std::vector<std::string> ret;
        ret.push_back( "n" );
        
        return ret;
    }

    static const std::vector<std::string> lowering_ops()
    {
        std::vector<std::string> ret;
        ret.push_back( "b" );
        
        return ret;
    }

    static const std::vector<std::string> raising_ops()
    {
        std::vector<std::string> ret;
        ret.push_back( "bdag" );
        
        return ret;
    }
};

struct TwoU1System {
    
    typedef TwoU1 grp;
    
    static const bool fermionic()   { return true; }
    static const int max_bond_dim() { return 20; }
    
    static const std::string lattice_lib() { return "alps";  }
    static const std::string model_lib()   { return "alps"; }
    
    static DmrgParameters parms()
    {
        DmrgParameters p;
        p.set("max_bond_dimension", max_bond_dim());
        p.set("lattice_library", lattice_lib());
        p.set("model_library", model_lib());
        p << model_parms();
        return p;
    }
    static ModelParameters model_parms()
    {
        ModelParameters p;
        p.set("LATTICE", "open chain lattice");
        p.set("L", 10);
        
        p.set("MODEL",                    "fermion Hubbard");
        p.set("Nmax",                     2                );
        p.set("t",                        1.               );
        p.set("U",                        1.               );
        p.set("CONSERVED_QUANTUMNUMBERS", "Nup,Ndown"      );
        p.set("Nup_total",                4                );
        p.set("Ndown_total",              4                );
        
        return p;
    }
    
    static const std::vector<std::string> dens_ops()
    {
        std::vector<std::string> ret;
        ret.push_back( "n_up" );
        ret.push_back( "n_down" );
        
        return ret;
    }
    
    static const std::vector<std::string> lowering_ops()
    {
        std::vector<std::string> ret;
        ret.push_back( "c_up" );
        ret.push_back( "c_down" );
        
        return ret;
    }
    
    static const std::vector<std::string> raising_ops()
    {
        std::vector<std::string> ret;
        ret.push_back( "cdag_up" );
        ret.push_back( "cdag_down" );
        
        return ret;
    }
};

struct U1DegSystem {
    
    typedef U1 grp;
    
    static const bool fermionic()   { return true; }
    static const int max_bond_dim() { return 20; }
    
    static const std::string lattice_lib() { return "alps";  }
    static const std::string model_lib()   { return "alps"; }
    
    static DmrgParameters parms()
    {
        DmrgParameters p;
        p.set("max_bond_dimension", max_bond_dim());
        p.set("lattice_library", lattice_lib());
        p.set("model_library", model_lib());
        p << model_parms();
        return p;
    }
    static ModelParameters model_parms()
    {
        ModelParameters p;
        p.set("LATTICE", "open chain lattice");
        p.set("L", 10);
        
        p.set("MODEL",           "alternative fermion Hubbard");
        p.set("t",               1.);
        p.set("U",               1.);
        p.set("CONSERVED_QUANTUMNUMBERS", "N");
        p.set("N_total",                  6);
        
        return p;
    }
    
    static const std::vector<std::string> dens_ops()
    {
        std::vector<std::string> ret;
        ret.push_back( "n_up" );
        ret.push_back( "n_down" );
        
        return ret;
    }
    
    static const std::vector<std::string> lowering_ops()
    {
        std::vector<std::string> ret;
        ret.push_back( "c_up" );
        ret.push_back( "c_down" );
        
        return ret;
    }
    
    static const std::vector<std::string> raising_ops()
    {
        std::vector<std::string> ret;
        ret.push_back( "cdag_up" );
        ret.push_back( "cdag_down" );
        
        return ret;
    }
};


template <class SymmGroup>
std::vector<double> measure_local(MPS<matrix, SymmGroup> const& mps,
                                  block_matrix<matrix, SymmGroup> const& ident,
                                  block_matrix<matrix, SymmGroup> const& op)
{
    typedef std::vector<block_matrix<matrix, SymmGroup> > op_vec;
    op_vec id_vec(1, ident);
    std::vector<double> vals(mps.size());
    boost::shared_ptr<lattice_impl> lat_ptr(new ChainLattice(mps.length()));
    Lattice lattice(lat_ptr);
//    std::cout << "ident is:\n" << ident << std::endl;
    for (int p=0; p<mps.size(); ++p) {
        generate_mpo::MPOMaker<matrix, SymmGroup> mpom(lattice, id_vec, id_vec);
        generate_mpo::Operator_Term<matrix, SymmGroup> term;
        term.operators.push_back( std::make_pair(p, op) );
        term.fill_operator = ident;
        mpom.add_term(term);
        MPO<matrix, SymmGroup> mpo = mpom.create_mpo();
        
        vals[p] = expval(mps, mpo);
    }
    return vals;
}

template <class SymmGroup>
double measure_ops(MPS<matrix, SymmGroup> const& mps,
                                block_matrix<matrix, SymmGroup> const& ident, block_matrix<matrix, SymmGroup> const& fill,
                                block_matrix<matrix, SymmGroup> const& op1, size_t pos1,
                                block_matrix<matrix, SymmGroup> const& op2, size_t pos2)
{
    typedef std::vector<block_matrix<matrix, SymmGroup> > op_vec;
    op_vec id_vec(1, ident);
    boost::shared_ptr<lattice_impl> lat_ptr(new ChainLattice(mps.length()));
    Lattice lattice(lat_ptr);
    generate_mpo::MPOMaker<matrix, SymmGroup> mpom(lattice, id_vec, id_vec);
    generate_mpo::Operator_Term<matrix, SymmGroup> term;
    block_matrix<matrix, SymmGroup> tmp;
    gemm(fill, op1, tmp);

    term.operators.push_back( std::make_pair(pos1, tmp) );
    for (int p=pos1+1; p<pos2; ++p) {
        term.operators.push_back( std::make_pair(p, fill) );
    }

    term.operators.push_back( std::make_pair(pos2, op2) );
    term.fill_operator = ident;
    mpom.add_term(term);
    MPO<matrix, SymmGroup> mpo = mpom.create_mpo();
    
    return expval(mps, mpo);
}


typedef boost::mpl::list<U1System, TwoU1System, U1DegSystem> test_systems;


BOOST_AUTO_TEST_CASE_TEMPLATE( dens_meas, ML, test_systems )
{
    typedef typename ML::grp grp;
    typedef typename grp::charge charge;
    
    DmrgParameters parms = ML::parms();
    
    Lattice lat(parms);
    Model<matrix,grp> model(lat, parms);
    
    block_matrix<matrix, grp> ident = model.identity_matrix();
    block_matrix<matrix, grp> fill  = model.filling_matrix();

    Index<grp> phys = model.phys_dim();
    std::cout << "phys: " << phys << std::endl;
    charge initc = model.total_quantum_numbers( parms );
    default_mps_init<matrix, grp> initializer(parms, std::vector<Index<grp> >(1, phys), initc, std::vector<int>(lat.size(),0));

    MPS<matrix, grp> mps(lat.size(), initializer);
    
    
    for (int i=0; i<ML::dens_ops().size(); ++i) {
        block_matrix<matrix, grp> dens   = model.get_operator( ML::dens_ops()[i] );
        block_matrix<matrix, grp> raise  = model.get_operator( ML::raising_ops()[i] );
        block_matrix<matrix, grp> lower  = model.get_operator( ML::lowering_ops()[i] );
        
        std::vector<double> meas1 = measure_local(mps, ident, dens);
        
        for (int it=0; it<3; ++it) {
            int n = drand48() * lat.size();
            maquis::cout << "measuring at n=" << n << std::endl;
            MPS<matrix, grp> mps2(mps);
            
            mps2.apply(lower, n);
            mps2.apply(raise, n);
            
            double meas = overlap(mps, mps2);
            
            BOOST_CHECK_CLOSE(meas1[n], meas, 1e-6 );
        }
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE( obdm_meas, ML, test_systems )
{
    typedef typename ML::grp grp;
    typedef typename grp::charge charge;
    
    DmrgParameters parms = ML::parms();
    
    Lattice lat(parms);
    Model<matrix,grp> model(lat, parms);
    
    block_matrix<matrix, grp> ident = model.identity_matrix();
    block_matrix<matrix, grp> fill  = model.filling_matrix();
    
    Index<grp> phys = model.phys_dim();
    std::cout << "phys: " << phys << std::endl;
    charge initc = model.total_quantum_numbers( parms );
    default_mps_init<matrix, grp> initializer(parms, std::vector<Index<grp> >(1, phys), initc, std::vector<int>(lat.size(),0));
    MPS<matrix, grp> mps(lat.size(), initializer);
    
    for (int i=0; i<ML::dens_ops().size(); ++i) {
        block_matrix<matrix, grp> dens   = model.get_operator( ML::dens_ops()[i] );
        block_matrix<matrix, grp> raise  = model.get_operator( ML::raising_ops()[i] );
        block_matrix<matrix, grp> lower  = model.get_operator( ML::lowering_ops()[i] );
        
        std::vector<double> meas1 = measure_local(mps, ident, dens);
        
        for (int it=0; it<3; ++it) {
            int ni = drand48() * (lat.size()-1);
            int nj = ni;
            while (nj <= ni)
                nj = drand48() * lat.size();
            maquis::cout << "measuring at ni=" << ni << ", nj=" << nj << std::endl;

            MPS<matrix, grp> mps2(mps);
            
            /// canonical order!
            if ( ML::fermionic() )
                mps2.apply(fill, lower, nj);
            else
                mps2.apply(lower, nj);

            if ( ML::fermionic() )
                mps2.apply(fill, raise, ni);
            else
                mps2.apply(raise, ni);
            

            double meas = overlap(mps, mps2);
            
            BOOST_CHECK_CLOSE(measure_ops(mps, ident, fill, raise, ni, lower, nj), meas, 1e-6 );
        }
    }
}


