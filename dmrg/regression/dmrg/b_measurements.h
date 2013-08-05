#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include "dmrg/mp_tensors/mps_mpo_ops.h"

#include "b_adjacency.h"
#include "b_generate_mpo.h"
#include "b_hamiltonians.h"

#include "b_DmrgParameters.h"

template<class Matrix, class SymmGroup, class Archive>
void measure_2pt_correlation(MPS<Matrix, SymmGroup> & mps,
                             b_adj::Adjacency & adj,
                             block_matrix<Matrix, SymmGroup> const & identity,
                             block_matrix<Matrix, SymmGroup> const & fill,
                             Archive & ar,
                             std::vector<block_matrix<Matrix, SymmGroup> > & ops,
                             std::string base_path)
{
    std::vector<double> dc;
    std::vector<std::string> labels;
    for (size_t p = 0; p < adj.size()-1; ++p) {
        b_mpos::CorrMaker<Matrix, SymmGroup> dcorr(mps.length(), identity, fill,
                                                   ops, p);
        MPO<Matrix, SymmGroup> mpo = dcorr.create_mpo();
        
        std::vector<double> dct = maquis::real((multi_expval(mps, mpo));
        std::copy(dct.begin(), dct.end(), std::back_inserter(dc));
        
        std::vector<std::string> lbt = dcorr.label_strings();
        std::copy(lbt.begin(), lbt.end(), std::back_inserter(labels));
    }
    
    ar[base_path + std::string("/mean/value")] << dc;
    ar[base_path + std::string("/labels")] << labels;
}
template<class Matrix, class SymmGroup, class Archive>
void measure_4ptnn_correlation(MPS<Matrix, SymmGroup> & mps,
                               b_adj::Adjacency & adj,
                               block_matrix<Matrix, SymmGroup> const & identity,
                               block_matrix<Matrix, SymmGroup> const & fill,
                               Archive & ar,
                               std::vector<block_matrix<Matrix, SymmGroup> > & ops,
                               std::string base_path)
{
    std::vector<double> dc;
    std::vector<std::string> labels;
    for (size_t p = 0; p < adj.size()-3; ++p) {
        b_mpos::CorrMakerNN<Matrix, SymmGroup> dcorr(mps.length(), identity, fill,
                                                     ops, p);
        MPO<Matrix, SymmGroup> mpo = dcorr.create_mpo();
        
        std::vector<double> dct = maquis::real(multi_expval(mps, mpo));
        std::copy(dct.begin(), dct.end(), std::back_inserter(dc));
        
        std::vector<std::string> lbt = dcorr.label_strings();
        std::copy(lbt.begin(), lbt.end(), std::back_inserter(labels));
    }
    
    ar[base_path + std::string("/mean/value")] << dc;
    ar[base_path + std::string("/labels")] << labels;
}

template<class Matrix, class SymmGroup>
struct measure_
{
    template<class Archive>
    void operator()(MPS<Matrix, SymmGroup> & mps,
                    b_adj::Adjacency & adj,
                    b_mpos::Hamiltonian<Matrix, SymmGroup> & H,
                    BaseParameters & model,
                    Archive & ar)
    { }
};

// these include specializations of measure_

#include "b_u1_measurements.h"
#include "b_2u1_measurements.h"

template<class Matrix, class SymmGroup, class Archive>
void measure(MPS<Matrix, SymmGroup> & mps,
             b_adj::Adjacency & adj,
             b_mpos::Hamiltonian<Matrix, SymmGroup> & H,
             BaseParameters & model,
             Archive & ar)
{
    measure_<Matrix, SymmGroup>()(mps, adj, H, model, ar);
}

#endif
