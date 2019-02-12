/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2018         Leon Freitag <lefreita@ethz.ch>
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
#include <fstream>
#include "dmrg/sim/matrix_types.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/mp_tensors/multi_canonize.h"
#include <boost/algorithm/string/split.hpp>
#include "dmrg/version.h"

/*
    Here we implement a simultaneous canonicalisation of MPS that have been
    optimised simultaneously
*/

#if defined(USE_TWOU1)
typedef TwoU1 grp;
#elif defined(USE_TWOU1PG)
typedef TwoU1PG grp;
#elif defined(USE_SU2U1)
typedef SU2U1 grp;
#elif defined(USE_SU2U1PG)
typedef SU2U1PG grp;
#elif defined(USE_NONE)
typedef TrivialGroup grp;
#elif defined(USE_U1)
typedef U1 grp;
#endif

int main(int argc, char ** argv)
{
    try {
        if (argc != 2)
        {
            maquis::cout << "Usage: <parms> <model_parms>" << std::endl;
            exit(1);
        }

        maquis::cout.precision(10);

        /// Loading parameters
        std::ifstream param_file(argv[1]);
        if (!param_file) {
            maquis::cerr << "Could not open parameter file." << std::endl;
            exit(1);
        }
        DmrgParameters parms(param_file);

        std::vector<std::string> checkpoints;
        std::string parm_checkpoint = parms["checkpoints"];
        boost::split(checkpoints, parm_checkpoint, boost::is_any_of(";"));

        std::vector<std::string> out_checkpoints;
        std::string parm_out_checkpoint = parms["out_checkpoints"];
        boost::split(out_checkpoints, parm_out_checkpoint, boost::is_any_of(";"));

        int site = parms["lrparam_site"];

        if(checkpoints.size() != out_checkpoints.size())
            throw std::runtime_error("Fatal error: number of output checkpoints must be the same as of input checkpoints");

        MPS<matrix, grp> mps;

        std::vector<MPS<matrix, grp> > mps_vec;

        for (auto&& chkp: checkpoints)
        {
            load(chkp, mps);
            // maquis::cout << "== Input MPS ==" << std::endl;
            // for (int i = 0; i < mps.length(); i++)
            //     maquis::cout << "MPS[" << i << "]" << std::endl << mps[i] << std::endl;

            mps_vec.emplace_back(std::move(mps));
        }

        multi_canonize(mps_vec, site);

        for (int i = 0; i < mps_vec.size(); i++)
        {
            save(out_checkpoints[i], mps_vec[i]);
            storage::archive ar(out_checkpoints[i]+"/props.h5", "w");
            ar["/parameters"] << parms;
            ar["/version"] << DMRG_VERSION_STRING;
            ar["/status/sweep"] << 0;
        }

        // for (auto&& mps_i: mps_vec)
        // {
        //     maquis::cout << "== Output MPS ==" << std::endl;
        //     for (int i = 0; i < mps_i.length(); i++)
        //         maquis::cout << "MPS[" << i << "]" << std::endl << mps_i[i] << std::endl;
        // }


    } catch (std::exception & e) {
        maquis::cerr << "Exception thrown!" << std::endl;
        maquis::cerr << e.what() << std::endl;
        exit(1);
    }
}