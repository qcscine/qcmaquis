/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2012-2012 by Sebastian Keller <sebkelle@phys.ethz.ch>
 * 
 *
 *****************************************************************************/

#ifndef MPS_INIT_HF_HPP
#define MPS_INIT_HF_HPP

#include "dmrg/mp_tensors/compression.h"

template<class Matrix>
struct hf_mps_init : public mps_initializer<Matrix, TwoU1>
{
    hf_mps_init(BaseParameters& params_) : params(params_) {}

    typedef Lattice::pos_t pos_t;
    typedef std::size_t size_t;

    void operator()(MPS<Matrix, TwoU1> & mps,
                    std::size_t Mmax,
                    Index<TwoU1> const & phys,
                    TwoU1::charge right_end)
    {
        default_mps_init<Matrix, TwoU1> di;
        di.init_sectors(mps, 5, phys, right_end, true);

        std::vector<std::size_t> hf_init = params.template get<std::vector<std::size_t> >("hf_occ");
        if (hf_init.size() != mps.length())
            throw std::runtime_error("HF occupation vector length != MPS length\n");

        TwoU1::charge max_charge = TwoU1::IdentityCharge;
        for (int i = 0; i < mps.length(); ++i)
        {
            mps[i].multiply_by_scalar(0.);

            size_t site_charge_up, site_charge_down, sc_input = hf_init[i];
            if (sc_input > 4)
                throw std::runtime_error(
                    "The hf_occ format has been changed to: 1=empty, 2=down, 3=up, 4=updown\n (not cumulative anymore)\n"
                );
            if (sc_input == 1 && i == 0)
                maquis::cout 
                    << "WARNING: The hf_occ format has been changed to: 1=empty, 2=down, 3=up, 4=updown\n (not cumulative anymore)\n";

            switch(sc_input) {
                case 4:
                    site_charge_up = 1;
                    site_charge_down = 1;
                    break;
                case 3:
                    site_charge_up = 1;
                    site_charge_down = 0;
                    break;
                case 2:
                    site_charge_up = 0;
                    site_charge_down = 1;
                    break;
                case 1:
                    site_charge_up = 0;
                    site_charge_down = 0;
                    break;
            }

            max_charge[0] += site_charge_up;
            max_charge[1] += site_charge_down;
            // Set largest charge sector = all 1
            size_t max_pos = mps[i].data().left_basis().position(max_charge);
            Matrix mfirst = mps[i].data()[max_pos];
            size_t nrow = mfirst.num_rows();
            size_t ncol = mfirst.num_cols();
            mps[i].data()[max_pos] = Matrix(nrow, ncol, 1.);

            mps[i].multiply_by_scalar(1. / mps[i].scalar_norm());

            //maquis::cout << std::endl;
            //maquis::cout << "mps[" << i << "]:\n" << mps[i] << std::endl;
        }

        //mps = compression::l2r_compress(mps, Mmax, 1e-6); 

        //maquis::cout << "\nMPS AFTER COMPRESSION:\n";
        //for(int i = 0; i < mps.length(); ++i) {
        //    maquis::cout << "mps[" << i << "]:\n" << mps[i] << std::endl;
        //    maquis::cout << mps[i].scalar_norm() << std::endl;
        //}
    }

    BaseParameters params;
};

#endif
