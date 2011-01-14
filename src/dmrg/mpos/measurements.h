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

        block_matrix<Matrix, U1> sz;
        sz.insert_block(Matrix(1, 1, 1), 1, 1);
        sz.insert_block(Matrix(1, 1, 0), 0, 0);
        sz.insert_block(Matrix(1, 1, -1), -1, -1);
        
        for (std::size_t p = 0; p < adj.size(); ++p)
        {
            mpos::MPOMaker<Matrix, U1> mpom(adj, H);
            std::vector<std::pair<std::size_t, block_matrix<Matrix, U1> > > v;
            v.push_back( std::make_pair( p, sz ) );
            mpom.add_term(v);
            MPO<Matrix, U1> mpo = mpom.create_mpo();
            
            double val = expval(mps, mpo, 0);
            magns.push_back(val);
        }

        ar << alps::make_pvp("/spectrum/results/Magnetization/mean/value", magns);
    }
    
    void operator()(MPS<Matrix, U1> & mps,
                    Adjacency & adj,
                    mpos::Hamiltonian<Matrix, U1> & H,
                    BaseParameters & model,
                    alps::hdf5::oarchive & ar)
    {
        if (model.get<std::string>("model") == std::string("biquadratic"))
            measure_blbq(mps, adj, H, model, ar);
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
