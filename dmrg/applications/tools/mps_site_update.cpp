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
#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/version.h"

/*
    Update one site of the MPS with the contents from file
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
        std::string mpstensorfile = parms["lrsite_update"];
        // site to update
        int site = parms["lrparam_site"];

        MPS<matrix, grp> mps;
        load(checkpoint, mps);
        typedef typename matrix::value_type value_type;
        std::vector<value_type> aux_elements;
        // read and parse the file
        std::ifstream infile(mpstensorfile);
        if (infile)
            std::copy(std::istream_iterator<value_type>(infile), std::istream_iterator<value_type>(), std::back_inserter(aux_elements));
        else
            throw std::runtime_error("File " + mpstensorfile + " could not be opened!");

        // Replace the MPSTensor at site site with the contents from file
        // TODO: do some basic checks

        MPSTensor<matrix, grp> & mpst = mps[site];
        mpst.make_left_paired();
        assert(mpst.data().num_elements() == aux_elements.size());

        size_t fileidx = 0;
        for (size_t i = 0; i < mpst.data().n_blocks(); i++)
        for (size_t j = 0; j < mpst.data()[i].num_rows(); j++)
        for (size_t k = 0; k < mpst.data()[i].num_cols(); k++)
        {
            parallel::guard::serial guard;
            mpst.data()[i](j,k) = aux_elements[fileidx++];
        }

        mpst.make_left_paired();

        // save the updated MPS
        save(checkpoint, mps);
        storage::archive ar(checkpoint+"/props.h5", "w");
        ar["/parameters"] << parms;
        ar["/version"] << DMRG_VERSION_STRING;
        ar["/status/sweep"] << 0;

    } catch (std::exception & e) {
        maquis::cerr << "Exception thrown!" << std::endl;
        maquis::cerr << e.what() << std::endl;
        exit(1);
    }

}