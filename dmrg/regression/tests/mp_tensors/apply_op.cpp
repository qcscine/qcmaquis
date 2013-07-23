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

#include "dmrg/utils/DmrgParameters2.h"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_initializers.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/models/generate_mpo.hpp"

#include "dmrg/models/factory.h"


typedef alps::numeric::matrix<double> matrix;


struct U1System {
    
    typedef U1 grp;
    
    static const bool fermionic()   { return false; }
    static const int max_bond_dim() { return 20; }
    
    static const std::string lattice_lib() { return "alps";  }
    static const std::string model_lib()   { return "coded"; }
    
    static ModelParameters parms()
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
    
    static ModelParameters parms()
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
    static const std::string model_lib()   { return "coded"; }
    
    static ModelParameters parms()
    {
        ModelParameters p;
        p.set("LATTICE", "open chain lattice");
        p.set("L", 10);
        
        p.set("MODEL",           "fermion Hubbard");
        p.set("t",               1.);
        p.set("U",               1.);
        p.set("u1_total_charge", 6);
        
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
    std::vector<double> vals(mps.size());
    for (int p=0; p<mps.size(); ++p) {
        generate_mpo::MPOMaker<matrix, SymmGroup> mpom(mps.length(), ident);
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
    generate_mpo::MPOMaker<matrix, SymmGroup> mpom(mps.length(), ident);
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
    
    ModelParameters parms = ML::parms();
    
    Lattice_ptr lat;
    typename model_traits<matrix, grp>::model_ptr model;
    model_parser<matrix, grp>(ML::lattice_lib(), ML::model_lib(), parms, lat, model);
    
    Index<grp> phys = model->get_phys();
    charge initc = model->initc( parms );
    default_mps_init<matrix, grp> initializer;
    MPS<matrix, grp> mps(lat->size(), ML::max_bond_dim(), phys, initc, initializer);
    
    block_matrix<matrix, grp> ident = model->get_op("id");
    block_matrix<matrix, grp> fill  = model->get_op("fill");
    
    for (int i=0; i<ML::dens_ops().size(); ++i) {
        block_matrix<matrix, grp> dens   = model->get_op( ML::dens_ops()[i] );
        block_matrix<matrix, grp> raise  = model->get_op( ML::raising_ops()[i] );
        block_matrix<matrix, grp> lower  = model->get_op( ML::lowering_ops()[i] );
        
        std::vector<double> meas1 = measure_local(mps, ident, dens);
        
        for (int it=0; it<3; ++it) {
            int n = drand48() * lat->size();
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
    
    ModelParameters parms = ML::parms();
    
    Lattice_ptr lat;
    typename model_traits<matrix, grp>::model_ptr model;
    model_parser<matrix, grp>(ML::lattice_lib(), ML::model_lib(), parms, lat, model);
    
    Index<grp> phys = model->get_phys();
    charge initc = model->initc( parms );
    default_mps_init<matrix, grp> initializer;
    MPS<matrix, grp> mps(lat->size(), ML::max_bond_dim(), phys, initc, initializer);
    
    block_matrix<matrix, grp> ident = model->get_op("id");
    block_matrix<matrix, grp> fill  = model->get_op("fill");
    
    for (int i=0; i<ML::dens_ops().size(); ++i) {
        block_matrix<matrix, grp> dens   = model->get_op( ML::dens_ops()[i] );
        block_matrix<matrix, grp> raise  = model->get_op( ML::raising_ops()[i] );
        block_matrix<matrix, grp> lower  = model->get_op( ML::lowering_ops()[i] );
        
        std::vector<double> meas1 = measure_local(mps, ident, dens);
        
        for (int it=0; it<3; ++it) {
            int ni = drand48() * (lat->size()-1);
            int nj = ni;
            while (nj <= ni)
                nj = drand48() * lat->size();
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


