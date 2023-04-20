/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef CHECKS_H_
#define CHECKS_H_


#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/utils/guess_symmetry.h"
#include "dmrg/mp_tensors/mps.h"

namespace maquis {
namespace checks {

    namespace detail{
        // get symmetry of a checkpoint as a string
        inline
        std::string get_symmetry(const std::string & chkpfile)
        {
            storage::archive ar_in(chkpfile+"/props.h5");
            BaseParameters chkp_parms;
            ar_in["/parameters"] >> chkp_parms;

            if (chkp_parms.defined("symmetry"))
                return chkp_parms["symmetry"].str();
            else
                throw std::runtime_error("Symmetry not defined in checkpoint "+chkpfile);
        }
    }
    // check whether the checkpoint chkpfile has the same symmetry as is declared in parms
    // removed guess_alps_symmetry as it is irrelevant for QCMaquis
    // Why is parms not const&, this could be potentially dangerous
    inline bool same_symmetry(BaseParameters & parms, const std::string & chkpfile)
    {
        return parms["symmetry"].str() == detail::get_symmetry(chkpfile);
    }

    inline void symmetry_check(BaseParameters & parms, const std::string & chkpfile)
    {
        std::string chkp_sym = detail::get_symmetry(chkpfile);
        if (parms["symmetry"].str() != chkp_sym)
            throw std::runtime_error("The existing checkpoint file " + chkpfile +  " has the wrong symmetry group " + chkp_sym + " instead of " + parms["symmetry"].str() + "\n");
    }

    // check whether the checkpoint has symmetry that ends in pg
    inline bool has_pg(const std::string & chkpfile)
    {
        std::string chkp_sym = detail::get_symmetry(chkpfile);
        std::regex regex_pg("pg$");

        return std::regex_search(chkp_sym, regex_pg);
    }


    template <class Matrix, class SymmGroup>
    void right_end_check(std::string filename, MPS<Matrix, SymmGroup> const & mps, typename SymmGroup::charge right_end)
    {
        if (right_end != mps[mps.size()-1].col_dim()[0].first) {
            std::stringstream parm_sector; parm_sector << right_end;
            std::stringstream mps_sector; mps_sector << mps[mps.size()-1].col_dim()[0].first;
            throw std::runtime_error("The existing checkpoint file " + filename +
                                     " has a different target symmetry sector (" + mps_sector.str() +
                                     ") compared to the input (" + parm_sector.str() + ")\n");
        }
    }

    template <class Matrix, class SymmGroup>
    void right_end_check(MPS<Matrix, SymmGroup> const & mps, typename SymmGroup::charge right_end)
    {
        if (right_end != mps[mps.size()-1].col_dim()[0].first) {
            std::stringstream parm_sector; parm_sector << right_end;
            std::stringstream mps_sector; mps_sector << mps[mps.size()-1].col_dim()[0].first;
            throw std::runtime_error("The given MPS has a different target symmetry sector (" + mps_sector.str() +
                                     ") compared to the input (" + parm_sector.str() + ")\n");
        }
    }



    inline void orbital_order_check(BaseParameters & parms, std::string chkpfile)
    {
        storage::archive ar_in(chkpfile+"/props.h5");
        BaseParameters chkp_parms;
        ar_in["/parameters"] >> chkp_parms;

        std::string chkp_order, parm_order;
        if (chkp_parms.defined("orbital_order")) {
            chkp_order = chkp_parms["orbital_order"].str();
            if (parms.is_set("orbital_order")) {
                parm_order = parms["orbital_order"].str();
                if (chkp_order != parm_order)
                    // make sure we apply the orbital order previously used for this MPS
                    parms.set("orbital_order",chkp_order);
                    //throw std::runtime_error("The existing checkpoint file " + chkpfile +  " has the wrong orbital order " + chkp_order + "\n");
            }
            else{
                // make sure we apply the orbital order previously used for this MPS
                parms.set("orbital_order",chkp_order);
            }
        }
    }
}
}


#endif
