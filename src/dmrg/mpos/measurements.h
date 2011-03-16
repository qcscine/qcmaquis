/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include <alps/hdf5.hpp>

#include "mp_tensors/mps_mpo_ops.h"

#include "mpos/adjacency.h"
#include "mpos/generate_mpo.h"
#include "mpos/hamiltonians.h"

#include "utils/DmrgParameters.h"

template<class Matrix, class SymmGroup>
struct measure_
{
    void operator()(MPS<Matrix, SymmGroup> & mps,
                    adj::Adjacency & adj,
                    mpos::Hamiltonian<Matrix, SymmGroup> & H,
                    BaseParameters & model,
                    alps::hdf5::oarchive & ar)
    { }
};

template<class Matrix, class SymmGroup>
void measure_2pt_correlation(MPS<Matrix, SymmGroup> & mps,
                             adj::Adjacency & adj,
                             block_matrix<Matrix, SymmGroup> const & identity,
                             block_matrix<Matrix, SymmGroup> const & fill,
                             alps::hdf5::oarchive & ar,
                             std::vector<block_matrix<Matrix, SymmGroup> > & ops,
                             std::string base_path)
{
    std::vector<double> dc;
    std::vector<std::string> labels;
    for (size_t p = 0; p < adj.size()-1; ++p) {
        mpos::CorrMaker<Matrix, U1> dcorr(mps.length(), identity, fill,
                                          ops, p);
        MPO<Matrix, U1> mpo = dcorr.create_mpo();
        
        std::vector<double> dct = multi_expval(mps, mpo);
        std::copy(dct.begin(), dct.end(), std::back_inserter(dc));
        
        std::vector<std::string> lbt = dcorr.label_strings();
        std::copy(lbt.begin(), lbt.end(), std::back_inserter(labels));
    }
    
    //std::copy(labels.begin(), labels.end(), std::ostream_iterator<std::string>(cout, " ")); cout << endl;
    //ar << alps::make_pvp(base_path + std::string("/labels"), labels);
    ar << alps::make_pvp(base_path + std::string("/mean/value"), dc);
    ar << alps::make_pvp(base_path + std::string("/labels"), labels);
}

template<class Matrix>
struct measure_<Matrix, U1>
{
    void measure_blbq(MPS<Matrix, U1> & mps,
                      adj::Adjacency & adj,
                      mpos::Hamiltonian<Matrix, U1> & H,
                      BaseParameters & model,
                      alps::hdf5::oarchive & ar)
    {
        std::vector<double> magns;

        std::vector<block_matrix<Matrix, U1> > ops(4);
        std::vector<std::string> names;
        
        block_matrix<Matrix, U1> ident = H.get_free();
        
        ops[0].insert_block(Matrix(1, 1, 1), 1, 1);
        ops[0].insert_block(Matrix(1, 1, 0), 0, 0);
        ops[0].insert_block(Matrix(1, 1, -1), -1, -1);
        names.push_back("Magnetization");
        
        ops[1].insert_block(Matrix(1, 1, 1), -1, -1);
        names.push_back("ColorDensity1");
        ops[2].insert_block(Matrix(1, 1, 1), 0, 0);
        names.push_back("ColorDensity2");
        ops[3].insert_block(Matrix(1, 1, 1), 1, 1);
        names.push_back("ColorDensity3");
        
        for (int i = 0; i < 4; ++i) {
            magns.clear();
            
            for (std::size_t p = 0; p < adj.size(); ++p)
            {
                mpos::MPOMaker<Matrix, U1> mpom(adj, H);
                std::vector<std::pair<int, block_matrix<Matrix, U1> > > v;
                v.push_back( std::make_pair( p, ops[i] ) );
                mpom.add_term(v);
                MPO<Matrix, U1> mpo = mpom.create_mpo();
                
                double val = expval(mps, mpo, 0);
                magns.push_back(val);
            }
            
            std::string n = std::string("/spectrum/results/") + names[i] + std::string("/mean/value");
            ar << alps::make_pvp(n, magns);
        }
     
        for (int i = 0; i < 4; ++i)
        {
            std::vector<block_matrix<Matrix, U1> > corr;
            corr.push_back( ops[i] );
            corr.push_back( ops[i] );
//            mpos::CorrMaker<Matrix, U1> mcorr(mps.length(), ident, ident, corr);
//            MPO<Matrix, U1> mpo = mcorr.create_mpo();
//            
//            std::vector<double> mc_v = multi_expval(mps, mpo);
            std::string name = std::string("/spectrum/results/") + names[i] + std::string("Correlation");
//            ar << alps::make_pvp(name + std::string("/labels"), mcorr.label_strings());
//            ar << alps::make_pvp(name + std::string("/mean/value"), mc_v);
            
            measure_2pt_correlation(mps, adj, ident, ident, ar,
                                    corr, name);
        }
        
    }
    
    void measure_superf(MPS<Matrix, U1> & mps,
                        adj::Adjacency & adj,
                        mpos::Hamiltonian<Matrix, U1> & H,
                        BaseParameters & model,
                        alps::hdf5::oarchive & ar)
    {
        std::vector<double> magns;
        
        block_matrix<Matrix, U1> dens;
        
        dens.insert_block(Matrix(1, 1, 1), 1, 1);
            
        for (std::size_t p = 0; p < adj.size(); ++p)
        {
            mpos::MPOMaker<Matrix, U1> mpom(adj, H);
            std::vector<std::pair<std::size_t, block_matrix<Matrix, U1> > > v;
            v.push_back( std::make_pair( p, dens ) );
            mpom.add_term(v);
            MPO<Matrix, U1> mpo = mpom.create_mpo();
            
            double val = expval(mps, mpo, 0);
            magns.push_back(val);
        }
            
        ar << alps::make_pvp("/spectrum/results/Density/mean/value", magns);
        
        std::vector<double> corrs;
        
        for (std::size_t p = 0; p < adj.size(); ++p)
        {
            std::vector<int> neighs = adj.forward(p);
            for (std::vector<int>::iterator it = neighs.begin();
                 it != neighs.end(); ++it)
            {
                mpos::MPOMaker<Matrix, U1> mpom(adj, H);
                std::vector<std::pair<std::size_t, block_matrix<Matrix, U1> > > v;
                v.push_back( std::make_pair( p, dens ) );
                v.push_back( std::make_pair(*it, dens) );
                mpom.add_term(v);
                MPO<Matrix, U1> mpo = mpom.create_mpo();
                
                double val = expval(mps, mpo, 0);
                corrs.push_back(val);
            }
        }
        
        ar << alps::make_pvp("/spectrum/results/NNDensityCorrelation/mean/value", corrs);
    }
    
    void measure_ff(MPS<Matrix, U1> & mps,
                    adj::Adjacency & adj,
                    mpos::Hamiltonian<Matrix, U1> & H,
                    BaseParameters & model,
                    alps::hdf5::oarchive & ar)
    {
        block_matrix<Matrix, U1> dens, create, destroy, sign, ident;
        
        dens.insert_block(Matrix(1, 1, 1), 1, 1);
        create.insert_block(Matrix(1, 1, 1), 0, 1);
        destroy.insert_block(Matrix(1, 1, 1), 1, 0);
        
        sign.insert_block(Matrix(1, 1, 1), 0, 0);
        sign.insert_block(Matrix(1, 1, -1), 1, 1);
        
        ident.insert_block(Matrix(1, 1, 1), 0, 0);
        ident.insert_block(Matrix(1, 1, 1), 1, 1);
        
        std::vector<double> density;
        for (std::size_t p = 0; p < adj.size(); ++p)
        {
            mpos::MPOMaker<Matrix, U1> mpom(adj, H);
            std::vector<std::pair<int, block_matrix<Matrix, U1> > > v;
            v.push_back( std::make_pair( p, dens ) );
            mpom.add_term(v);
            MPO<Matrix, U1> mpo = mpom.create_mpo();
            
            double val = expval(mps, mpo, 0);
            density.push_back(val);
        }
        ar << alps::make_pvp("/spectrum/results/Density/mean/value", density);
        
        std::vector<block_matrix<Matrix, U1> > density_corr;
        density_corr.push_back( dens );
        density_corr.push_back( dens );
        measure_2pt_correlation(mps, adj, ident, ident,
                                ar, density_corr,
                                "/spectrum/results/DensityCorrelation");
        
        std::vector<block_matrix<Matrix, U1> > onebody;
        onebody.push_back( create );
        onebody.push_back( destroy );
        measure_2pt_correlation(mps, adj, ident, sign,
                                ar, onebody,
                                "/spectrum/results/OneBodyDM");
    }
    
    
    void operator()(MPS<Matrix, U1> & mps,
                    adj::Adjacency & adj,
                    mpos::Hamiltonian<Matrix, U1> & H,
                    BaseParameters & model,
                    alps::hdf5::oarchive & ar)
    {
        if (model.get<std::string>("model") == std::string("biquadratic"))
            measure_blbq(mps, adj, H, model, ar);
        else if (model.get<std::string>("model") == std::string("FreeFermions"))
            measure_ff(mps, adj, H, model, ar);
    }
};

template<class Matrix, class SymmGroup>
void measure(MPS<Matrix, SymmGroup> & mps,
             adj::Adjacency & adj,
             mpos::Hamiltonian<Matrix, SymmGroup> & H,
             BaseParameters & model,
             alps::hdf5::oarchive & ar)
{
    measure_<Matrix, SymmGroup>()(mps, adj, H, model, ar);
}

#endif
