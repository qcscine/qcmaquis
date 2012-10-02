//
//  smrg_strcorr.h
//  MAQUIS_DMRG
//
//  Created by Michele Dolfi on 7/23/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#ifndef MAQUIS_DMRG_STR_CORR_H
#define MAQUIS_DMRG_STR_CORR_H

#include "dmrg/block_matrix/detail/alps.hpp"
#ifdef USE_COMPLEX
typedef alps::numeric::matrix<std::complex<double> > matrix;
#else
typedef alps::numeric::matrix<double> matrix;
#endif

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

#include "dmrg/utils/DmrgParameters2.h"

#include <alps/hdf5.hpp>
#include <alps/numeric/real.hpp>
#include <boost/filesystem.hpp>


class StrCorr {
public:
    typedef block_matrix<matrix, U1> op_t;
    
    StrCorr(DmrgParameters & parms_, ModelParameters & model_)
    : parms(parms_)
    , model(model_)
    {
        // Loading state
        if (! boost::filesystem::exists(parms.get<std::string>("chkpfile")))
            throw std::runtime_error("Checkpoint file doesn not exist.");
        size_t graining = 0;
        {
            alps::hdf5::archive ar(parms.get<std::string>("chkpfile"));
            ar >> alps::make_pvp("/state", mps);
            if (ar.is_data("/status/graining"))
                ar >> alps::make_pvp("/status/graining", graining);
        }
        
        // Init model
        model = model_.get_at_index("graining", graining);
        init_model();
        L = model.get<int>("L"); N = model.get<int>("Ndiscr");
        Ltot = L * N;
        
        if (Ltot != mps.length())
            throw std::runtime_error("MPS doesn't match the continuum model.");
    }
    
    // Single-site string operator of length l grid points (from middle)
    // hardcoding unit filling
    void measure_ss_string_unit(size_t l)
    {
        op_t strop;
        strop.insert_block(matrix(1, 1, -1), 0, 0);
        for (int n=1; n<=model.get<int>("Nmax"); ++n)
        {
            if ((n-1) % 2 == 0)
                strop.insert_block(matrix(1, 1, 1), n, n);
            else
                strop.insert_block(matrix(1, 1, -1), n, n);
        }
        
        // build MPO
        MPO<matrix, U1> mpo(Ltot);
        size_t middle = Ltot / 2;
        for (size_t i=0; i<Ltot; ++i) {
            if (i >= middle - l/2 && i < middle + l/2)
                mpo[i](0,0) = strop;
            else
                mpo[i](0,0) = ident;
        }
        
        // eval & save
        matrix::value_type val = expval(mps,mpo);
        {
            alps::hdf5::archive ar(parms.get<std::string>("resultfile"), alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            ar << alps::make_pvp(std::string("/spectrum/results/String order SS ") + boost::lexical_cast<std::string>(l) + std::string("/mean/value"),
                                 std::vector<double>(1, maquis::real(val)));
        }
    }
    
    // Unite-cell string operator of length l unit cells
    void measure_uc_string(size_t l)
    {
        int filling = model.get<int>("u1_total_charge") / L;
        op_t strop;
        strop.insert_block(matrix(1, 1, 1), 0, 0);
        for (int n=1; n<=model.get<int>("Nmax"); ++n)
        {
            if (n % 2 == 0)
                strop.insert_block(matrix(1, 1, 1), n, n);
            else
                strop.insert_block(matrix(1, 1, -1), n, n);
        }
        
        // build MPO
        MPO<matrix, U1> mpo(Ltot);
        size_t middle = Ltot / 2;
        for (size_t i=0; i<Ltot; ++i) {
            if (i >= middle - l/2 && i < middle + l/2)
                mpo[i](0,0) = strop;
            else 
                mpo[i](0,0) = ident;
        }
        
        // eval & save
        matrix::value_type val = expval(mps,mpo) * exp(-filling*double(l/N));
        {
            alps::hdf5::archive ar(parms.get<std::string>("resultfile"), alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            ar << alps::make_pvp(std::string("/spectrum/results/String order UC ") + boost::lexical_cast<std::string>(l) + std::string("/mean/value"),
                                 std::vector<double>(1, maquis::real(val)));
        }
    }
    
    
private:
    
    void init_model()
    {
        phys.insert(std::make_pair(0, 1));
        ident.insert_block(matrix(1, 1, 1), 0, 0);
        
        for (int n=1; n<=model.get<int>("Nmax"); ++n)
        {
            phys.insert(std::make_pair(n, 1));
            ident.insert_block(matrix(1, 1, 1), n, n);
            count.insert_block(matrix(1, 1, n), n, n);
        }
    }
    
    
    DmrgParameters & parms;
    ModelParameters & model;
    size_t L, N, Ltot;
    
    Index<U1> phys;
    op_t ident, count;
    MPS<matrix, U1> mps;
};


#endif
