/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2015-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef CHECKS_HPP_
#define CHECKS_HPP_


#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/mp_tensors/mps.h"

std::string guess_alps_symmetry(BaseParameters & parms);

namespace maquis {
namespace checks {

    inline void symmetry_check(BaseParameters & parms, std::string chkpfile)
    {
        storage::archive ar_in(chkpfile+"/props.h5");
        BaseParameters chkp_parms;
        ar_in["/parameters"] >> chkp_parms;

        std::string chkp_sym, parm_sym;
        if (chkp_parms.defined("symmetry")) {
            chkp_sym = chkp_parms["symmetry"].str();
            parm_sym = parms["symmetry"].str();
        }
        else {
            chkp_sym = guess_alps_symmetry(chkp_parms);
            parm_sym = guess_alps_symmetry(parms);
        }

        if (chkp_sym != parm_sym)
            throw std::runtime_error("The existing checkpoint file " + chkpfile +  " has the wrong symmetry group " + chkp_sym + "\n");
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

}
}


#endif
