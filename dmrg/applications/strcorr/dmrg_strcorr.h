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

#include "dmrg/utils/DmrgParameters.h"

#include <alps/numeric/real.hpp>
#include <boost/filesystem.hpp>


class StrCorr {
public:
    typedef block_matrix<matrix, U1> op_t;
    
    StrCorr(DmrgParameters & parms_)
    : parms(parms_)
    {
        // Loading state
        if (! boost::filesystem::exists(parms["chkpfile"].str()))
            throw std::runtime_error("Checkpoint file doesn not exist.");
        size_t graining = 0;
        load(parms["chkpfile"].str(), mps);
        {
            storage::archive ar(parms["chkpfile"].str()+"/props.h5");
            if (ar.is_data("/status/graining"))
                ar["/status/graining"] >> graining;
        }
        
        // Init model
        parms = parms_.get_at_index("graining", graining);
        init_model();
        L = parms["L"]; N = parms["Ndiscr"];
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
        for (int n=1; n<=parms["Nmax"]; ++n)
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
                mpo[i].set(0,0, strop);
            else
                mpo[i].set(0,0, ident);
        }
        
        // eval & save
        matrix::value_type val = expval(mps,mpo);
        {
            storage::archive ar(parms["resultfile"].str(), "w");
            ar[std::string("/spectrum/results/String order SS ") + boost::lexical_cast<std::string>(l) + std::string("/mean/value")] << std::vector<double>(1, maquis::real(val));
        }
    }
    
    // Unite-cell string operator of length l unit cells
    void measure_uc_string(size_t l)
    {
        int filling = parms["u1_total_charge"] / L;
        op_t strop;
        strop.insert_block(matrix(1, 1, 1), 0, 0);
        for (int n=1; n<=parms["Nmax"]; ++n)
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
                mpo[i].set(0,0, strop);
            else 
                mpo[i].set(0,0, ident);
        }
        
        // eval & save
        matrix::value_type val = expval(mps,mpo) * exp(-filling*double(l/N));
        {
            storage::archive ar(parms["resultfile"].str(), "w");
            ar[std::string("/spectrum/results/String order UC ") + boost::lexical_cast<std::string>(l) + std::string("/mean/value")] << std::vector<double>(1, maquis::real(val));
        }
    }
    
    
private:
    
    void init_model()
    {
        phys.insert(std::make_pair(0, 1));
        ident.insert_block(matrix(1, 1, 1), 0, 0);
        
        for (int n=1; n<=parms["Nmax"]; ++n)
        {
            phys.insert(std::make_pair(n, 1));
            ident.insert_block(matrix(1, 1, 1), n, n);
            count.insert_block(matrix(1, 1, n), n, n);
        }
    }
    
    
    DmrgParameters & parms;
    size_t L, N, Ltot;
    
    Index<U1> phys;
    op_t ident, count;
    MPS<matrix, U1> mps;
};


#endif
