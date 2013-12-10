/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 * 
 *
 *****************************************************************************/

#ifndef MPS_INIT_HF_HPP
#define MPS_INIT_HF_HPP

#include "dmrg/mp_tensors/compression.h"
#include "dmrg/models/chem/pg_util.h"

template<class Matrix, class SymmGroup>
struct hf_mps_init : public mps_initializer<Matrix, SymmGroup>
{
    hf_mps_init(BaseParameters parms_, BaseParameters model_,
                std::vector<Index<SymmGroup> > const& phys_dims,
                typename SymmGroup::charge right_end,
                std::vector<int> const& site_type)
    : parms(parms_), model(model_)
    , di(parms, model, phys_dims, right_end, site_type)
    {}

    typedef Lattice::pos_t pos_t;
    typedef std::size_t size_t;

    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
        di.init_sectors(mps, 5, true, 0);

        std::vector<std::size_t> hf_init = model["hf_occ"];

        std::vector<pos_t> order(mps.length());
        if (!model.is_set("orbital_order"))
            for (pos_t p = 0; p < mps.length(); ++p)
                order[p] = p;
        else
            order = model["orbital_order"].template as<std::vector<pos_t> >();

        std::transform(order.begin(), order.end(), order.begin(), boost::lambda::_1-1);

        if (hf_init.size() != mps.length())
            throw std::runtime_error("HF occupation vector length != MPS length\n");

        std::vector<typename SymmGroup::subcharge> irreps = parse_symm<SymmGroup>(mps.length(), model);

        typename SymmGroup::charge max_charge = SymmGroup::IdentityCharge;
        for (pos_t i = 0; i < mps.length(); ++i)
        {
            mps[i].multiply_by_scalar(0.0);

            size_t sc_input = hf_init[order[i]];
            typename SymmGroup::charge site_charge(0);

            if (sc_input > 4)
                throw std::runtime_error(
                    "The hf_occ format has been changed to: 1=empty, 2=down, 3=up, 4=updown\n (not cumulative anymore)\n"
                );

            switch(sc_input) {
                case 4:
                    site_charge[0] = 1; // up
                    site_charge[1] = 1; // down
                    break;
                case 3:
                    site_charge[0] = 1;
                    site_charge[1] = 0;
                    PGCharge<SymmGroup>()(site_charge, irreps[i]);
                    break;
                case 2:
                    site_charge[0] = 0;
                    site_charge[1] = 1;
                    PGCharge<SymmGroup>()(site_charge, irreps[i]);
                    break;
                case 1:
                    site_charge[0] = 0;
                    site_charge[1] = 0;
                    break;
            }

            max_charge = SymmGroup::fuse(max_charge, site_charge);

            // Set largest charge sector = all 1
            size_t max_pos = mps[i].data().left_basis().position(max_charge);
            Matrix & mfirst = mps[i].data()[max_pos];
            size_t nrow = mfirst.num_rows();
            size_t ncol = mfirst.num_cols();
            mps[i].data()[max_pos] = Matrix(nrow, ncol, 1.);

            mps[i].multiply_by_scalar(1. / mps[i].scalar_norm());

        }

        //mps = compression::l2r_compress(mps, Mmax, 1e-6); 

        //maquis::cout << "\nMPS AFTER COMPRESSION:\n";
        //for(int i = 0; i < mps.length(); ++i) {
        //    maquis::cout << "mps[" << i << "]:\n" << mps[i] << std::endl;
        //    maquis::cout << mps[i].scalar_norm() << std::endl;
        //}
    }

    BaseParameters parms;
    BaseParameters model;
    
    default_mps_init<Matrix, SymmGroup> di;
};

#endif
