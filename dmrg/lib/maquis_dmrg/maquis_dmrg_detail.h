/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2018-2019    Leon Freitag <lefreita@ethz.ch>
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
#ifndef MAQUIS_DMRG_DETAIL_H
#define MAQUIS_DMRG_DETAIL_H

#include "integral_interface.h"
#include "dmrg/utils/results_collector.h"
#include "dmrg/utils/DmrgParameters.h"

namespace maquis {
    using chem::integral_map;
    using chem::serialize;
    // Constants for integral_map template specialisation, whether one should use relativistic or nonrelativistic integrals
    namespace integrals {
        const bool relativistic = true;
        const bool nonrelativistic = false;
    }

    namespace interface_detail {
        // Generate names for SU2U1 checkpoint files
        std::string su2u1_name(const std::string & pname, int state);

        // Generate names for 2U1 checkpoint files
        // Generates only the checkpoint with minimum Ms, as consistently used in the interface
        // returns a tuple containing the filename, Nup and Ndown electrons
        // If Ms is set to 0, the Ms closest to 0 (i.e. 0 or 1) will be chosen
        // (to allow for default values where we do not care for Ms)
        std::tuple<std::string, int, int> twou1_name_Nup_Ndown(const std::string & pname, int state, int nel, int multiplicity, int Ms=0);

        // same as above, but returns only the filename
        std::string twou1_name(const std::string & pname, int state, int nel, int multiplicity, int Ms=0);

        // SU2U1 result name file
        // TODO: Similar functions for (SU2U1 and 2U1) checkpoint names are found in mpssi_interface.cpp, maybe they should be moved here?
        inline std::string su2u1_result_name(const std::string& pname, int state)
        {
            std::string ret = pname + ".results_state." + std::to_string(state) + ".h5";
            return ret;
        }

        // Result file name for transition 3-RDM measurements
        inline std::string trans3rdm_result_name(const std::string& pname, int state, int bra_state)
        {
            std::string ret = pname + ".trans3rdm." + std::to_string(state) + "_" + std::to_string(bra_state) + ".h5";
            return ret;
        }
    }

    // Set parameters required for relativistic calculation
    void prepare_relativistic(BaseParameters& parms, bool magnetic = false);

    // Transforms SU2 checkpoint to 2U1 checkpoint
    void transform(const std::string & pname, int state, int Ms=0);

}


#endif