#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include <alps/hdf5.hpp>

#include "mp_tensors/mps_mpo_ops.h"

#include "mpos/adjancency.h"
#include "mpos/generate_mpo.h"
#include "mpos/hamiltonians.h"

#include "utils/DmrgParameters.h"

template<class Matrix, class SymmGroup>
struct measure_
{
    void operator()(MPS<Matrix, SymmGroup> & mps,
                    Adjacency & adj,
                    mpos::Hamiltonian<Matrix, SymmGroup> & H,
                    BaseParameters & model,
                    alps::hdf5::oarchive & ar)
    { }
};

template<class Matrix>
struct measure_<Matrix, U1>
{
    void measure_blbq(MPS<Matrix, U1> & mps,
                      Adjacency & adj,
                      mpos::Hamiltonian<Matrix, U1> & H,
                      BaseParameters & model,
                      alps::hdf5::oarchive & ar)
    {
        std::vector<double> magns;

        std::vector<block_matrix<Matrix, U1> > ops(4);
        std::vector<std::string> names;
        
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
                std::vector<std::pair<std::size_t, block_matrix<Matrix, U1> > > v;
                v.push_back( std::make_pair( p, ops[i] ) );
                mpom.add_term(v);
                MPO<Matrix, U1> mpo = mpom.create_mpo();
                
                double val = expval(mps, mpo, 0);
                magns.push_back(val);
            }
            
            std::string n = std::string("/spectrum/results/") + names[i] + std::string("/mean/value");
            ar << alps::make_pvp(n, magns);
        }
    }
    
    void measure_ff(MPS<Matrix, U1> & mps,
                    Adjacency & adj,
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
        
        std::vector<block_matrix<Matrix, U1> > density_corr;
        density_corr.push_back( dens );
        density_corr.push_back( dens );
        mpos::CorrMaker<Matrix, U1> dcorr(mps.length(), ident, ident,
                                                 density_corr);
        MPO<Matrix, U1> mpo = dcorr.create_mpo();
        
        std::vector<double> dc = multi_expval(mps, mpo);
        ar << alps::make_pvp("/spectrum/results/DensityCorrelation/labels", dcorr.label_strings());
        ar << alps::make_pvp("/spectrum/results/DensityCorrelation/mean/value", dc);
        
        std::vector<block_matrix<Matrix, U1> > onebody;
        onebody.push_back( create );
        onebody.push_back( destroy );
        mpos::CorrMaker<Matrix, U1> ob(mps.length(), ident, sign,
                                       onebody);
        MPO<Matrix, U1> mpo2 = ob.create_mpo();
        
        dc = multi_expval(mps, mpo2);
        ar << alps::make_pvp("/spectrum/results/OneBodyDM/labels", ob.label_strings());
        ar << alps::make_pvp("/spectrum/results/OneBodyDM/mean/value", dc);
    }
    
    
    void operator()(MPS<Matrix, U1> & mps,
                    Adjacency & adj,
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
             Adjacency & adj,
             mpos::Hamiltonian<Matrix, SymmGroup> & H,
             BaseParameters & model,
             alps::hdf5::oarchive & ar)
{
    measure_<Matrix, SymmGroup>()(mps, adj, H, model, ar);
}

#endif
