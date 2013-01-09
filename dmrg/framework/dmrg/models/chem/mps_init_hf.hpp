/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2012-2012 by Sebastian Keller <sebkelle@phys.ethz.ch>
 * 
 *
 *****************************************************************************/

#ifndef MPS_INIT_HF_HPP
#define MPD_INIT_HF_HPP

#include "dmrg/mp_tensors/compression.h"

template<class Matrix>
struct hf_mps_init : public mps_initializer<Matrix, TwoU1>
{
    hf_mps_init(BaseParameters& params_) : params(params_) {}

    void operator()(MPS<Matrix, TwoU1> & mps,
                    std::size_t Mmax,
                    Index<TwoU1> const & phys,
                    TwoU1::charge right_end)
    {
        default_mps_init<Matrix, TwoU1> di;
        di.init_sectors(mps, 5, phys, right_end, true);

        for (int i = 0; i < mps.length(); ++i)
        {
            mps[i].multiply_by_scalar(0.);

            std::vector<std::size_t> hf_init = params.template get<std::vector<std::size_t> >("hf_occ");
            if (hf_init.size() != mps.length()) {
                std::cout << "HF occupation vector length != MPS length\n";
                exit(1);
            }

            TwoU1::charge max_charge;
            //max_charge[0] = std::min(i+1, 4);
            max_charge[0] = hf_init[i];
            max_charge[1] = hf_init[i];
            // Set largest charge sector = all 1
            std::size_t max_pos = mps[i].data().left_basis().position(max_charge);
            Matrix mfirst = mps[i].data()[max_pos];
            std::size_t nrow = mfirst.num_rows();
            std::size_t ncol = mfirst.num_cols();
            mps[i].data()[max_pos] = Matrix(nrow, ncol, 1.);

            mps[i].multiply_by_scalar(1. / mps[i].scalar_norm());

            //std::cout << std::endl;
            //std::cout << "mps[" << i << "]:\n" << mps[i] << std::endl;
        }

        //mps = compression::l2r_compress(mps, Mmax, 1e-6); 

        //std::cout << "\nMPS AFTER COMPRESSION:\n";
        //for(int i = 0; i < mps.length(); ++i) {
        //    std::cout << "mps[" << i << "]:\n" << mps[i] << std::endl;
        //    std::cout << mps[i].scalar_norm() << std::endl;
        //}
    }

    BaseParameters params;
};

#endif
