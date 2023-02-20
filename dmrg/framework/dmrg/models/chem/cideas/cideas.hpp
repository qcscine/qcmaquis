/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
#ifndef CIDEAS_HPP
#define CIDEAS_HPP

#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/models/chem/mps_init_deas.hpp"

namespace maquis
{
    // perform CI-DEAS
    template <class Matrix, class SymmGroup>
    MPS<Matrix, SymmGroup> cideas(DmrgParameters parms, const Matrix& s1)
    {
        Lattice lat(parms);
        Model<Matrix,SymmGroup> model(lat,parms);

        //create symmetry vector -> site_types
        std::vector<typename SymmGroup::subcharge> site_types(lat.size()), site_types_distinct;
        for (int i = 0; i < lat.size(); i++)
        {
            site_types[i] = lat.get_prop<typename SymmGroup::subcharge>("type", i);
            if(std::find(site_types_distinct.begin(),site_types_distinct.end(),site_types[i]) == site_types_distinct.end())
                site_types_distinct.push_back(site_types[i]);
        }

        std::cout <<"site_types: ";
        for(int i = 0; i < site_types.size(); i++)
            std::cout << site_types[i] <<", ";

        std::cout << std::endl;

        std::vector<Index<SymmGroup> > phys_dims;
        for (int i = 0; i < site_types_distinct.size(); i++)
            //phys_dims.push_back(model.phys_dim(site_types_distinct[i]));
            phys_dims.push_back(model.phys_dim(i));

        for (int i = 0; i < phys_dims.size(); i++)
            std::cout << "phys_dims["<<i<<"] = " << phys_dims[i] << std::endl;

        //get right end charge
        typename SymmGroup::charge right_end = chem::detail::qn_helper<SymmGroup>().total_qn(parms);
        maquis::cout << "Right end: " << right_end <<std::endl;

        //create MPS
        MPS<Matrix, SymmGroup> mps(lat.size());
        deas_mps_init<Matrix, SymmGroup> deas_creator(parms, s1, phys_dims, right_end, site_types);
        deas_creator(mps);

        return mps;
    }
}

#endif