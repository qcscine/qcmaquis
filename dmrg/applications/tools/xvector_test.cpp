/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2019         Leon Freitag <lefreita@ethz.ch>
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
#include "dmrg/mp_tensors/xvector.h"
#include "dmrg/utils/DmrgParameters.h"

/*
    Test for the transformation to non-redundant MPS parameters X
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

        std::string checkpoint = parms["chkpfile"];
        std::string checkpoint_out = parms["outchkpfile"];

        MPS<matrix, grp> mps;
        load(checkpoint, mps);

        lr::XVector<matrix, grp> x(mps, mps); // this should yield zeros

        maquis::cout << "Number of nonredundant MPS parameters: " << std::accumulate(x.begin(), x.end(), 0,
            [&](std::size_t sum, block_matrix<matrix, grp> m) {
                    for (int i = 0; i < m.n_blocks(); i++)
                        sum += m[i].num_cols()*m[i].num_rows();
                    return sum;
                }) << std::endl;

        x.save(checkpoint_out);
        if(parms.is_set("xvec_aux_file"))
            x.dump_to_textfile(parms["xvec_aux_file"]);

        lr::XVector<matrix, grp> x2;
        x2.load(checkpoint_out);
        if(parms.is_set("xvec_aux_file"))
            x2.update_from_textfile(parms["xvec_aux_file"]);

        // MPS<matrix, grp> mps_transformed = x.transformXtoB();
        // save(checkpoint_out, mps_transformed);
    } catch (std::exception & e) {
        maquis::cerr << "Exception thrown!" << std::endl;
        maquis::cerr << e.what() << std::endl;
        exit(1);
    }

}