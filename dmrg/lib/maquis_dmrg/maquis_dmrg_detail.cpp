/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2020         Leon Freitag <lefreita@ethz.ch>
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

#include "maquis_dmrg_detail.h"
#include "dmrg/models/chem/transform_symmetry.hpp"

namespace maquis {
    namespace interface_detail {

        // Generate names for SU2U1 checkpoint files
        std::string su2u1_name(const std::string & pname, int state)
        {
            // TODO: should we check if pname and state are allowed here too with allowed_names_states()?
            // if (allowed_names_states(pname, state) == -1)
            //     throw std::runtime_error("Do not know the project name" + pname );

            std::string ret = pname + ".checkpoint_state." + std::to_string(state) + ".h5";
            return ret;
        }

        // Generate names for 2U1 checkpoint files (helper)
        std::tuple<std::string, int, int>
        twou1_name_Nup_Ndown(const std::string & pname, int state, int nel, int multiplicity, int Ms)
        {
            int Nup = (nel + Ms)/2;
            int Ndown = (nel - Ms)/2;

            int remainder = Ms;
            if (Ms == 0)
            {
                // Use 2U1 checkpoint with Ms=0 or 1 for default Ms
                remainder = nel % 2; // integer division
                Nup = nel / 2 + remainder;
                Ndown = nel / 2 - remainder;
            }

            // Use 2U1 checkpoint with Ms=0 or 1

            // int remainder = nel % 2; // integer division
            // int Nup = nel / 2 + remainder;
            // int Ndown = nel / 2 - remainder;


            std::string ret = pname + ".checkpoint_state." + std::to_string(state)
                                    + "." + std::to_string(multiplicity) + "." + std::to_string(remainder)
                                    + ".h5";
            return std::make_tuple(ret, Nup, Ndown);
        }

        std::string twou1_name(const std::string & pname, int state, int nel, int multiplicity, int Ms)
        {
            return std::get<0>(twou1_name_Nup_Ndown(pname, state, nel, multiplicity, Ms));
        }

    }

// Transforms SU2 checkpoint to 2U1 checkpoint
// Mostly copy-paste from mps_transform.cpp, but creates only one 2U1 checkpoint per state
// corresponding to the state with the highest Sz

    void transform(const std::string & pname, int state, int Ms)
    {
        // This works only for double
        // Do we really need complex? (since we don't transform double groups)
        typedef alps::numeric::matrix<double> matrix;

    #if defined(HAVE_SU2U1PG)
                typedef SU2U1PG SU2U1grp;
                typedef TwoU1PG TwoU1grp;
    #elif defined(HAVE_SU2U1)
                typedef SU2U1 SU2U1grp;
                typedef TwoU1 TwoU1grp;
    #endif

        std::string checkpoint_name = interface_detail::su2u1_name(pname, state);

        BaseParameters parms;

        if (!boost::filesystem::exists(checkpoint_name))
            throw std::runtime_error("input MPS " + checkpoint_name + " does not exist\n");

        // load source MPS
        MPS<matrix, SU2U1grp> mps;
        load(checkpoint_name, mps);

        // fetch parameters and modify symmetry
        storage::archive ar_in(checkpoint_name + "/props.h5");

        ar_in["/parameters"] >> parms;
        parms.set("init_state", "const");
    #if defined(HAVE_SU2U1PG)
        parms.set("symmetry", "2u1pg");
    #elif defined(HAVE_SU2U1)
        parms.set("symmetry", "2u1");
    #endif
        int Nup, Ndown;
        std::string twou1_checkpoint_name;
        int nel = parms["nelec"];
        int multiplicity = parms["spin"];

        // get number of up/down electrons and the checkpoint name for the 2U1 checkpoint
        std::tie(twou1_checkpoint_name, Nup, Ndown) = maquis::interface_detail::twou1_name_Nup_Ndown(pname, state, nel, multiplicity, Ms);

        parms.set("u1_total_charge1", Nup);
        parms.set("u1_total_charge2", Ndown);

        // transform MPS
        MPS<matrix, TwoU1grp> mps_out = transform_mps<matrix, SU2U1grp>()(mps, Nup, Ndown);

        save(twou1_checkpoint_name, mps_out);

        if (boost::filesystem::exists(twou1_checkpoint_name + "/props.h5"))
            boost::filesystem::remove(twou1_checkpoint_name + "/props.h5");
        boost::filesystem::copy(checkpoint_name + "/props.h5", twou1_checkpoint_name + "/props.h5");

        storage::archive ar_out(twou1_checkpoint_name + "/props.h5", "w");
        ar_out["/parameters"] << parms;
    }

}