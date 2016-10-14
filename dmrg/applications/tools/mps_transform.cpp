/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Laboratory of Physical Chemistry, ETH Zurich
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

#ifdef USE_AMBIENT
#include <mpi.h>
#endif
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <sys/stat.h>

#include <boost/filesystem.hpp>

using std::cerr;
using std::cout;
using std::endl;

#ifdef USE_AMBIENT
    #include "dmrg/block_matrix/detail/ambient.hpp"
    typedef ambient::numeric::tiles<ambient::numeric::matrix<double> > Matrix;
#else
    #include "dmrg/block_matrix/detail/alps.hpp"
    typedef alps::numeric::matrix<double> Matrix;
#endif

#include "dmrg/models/model.h"
#include "dmrg/models/chem/transform_symmetry.hpp"

#if defined(USE_SU2U1)
typedef SU2U1 grp;
typedef TwoU1 mapgrp;
#elif defined(USE_SU2U1PG)
typedef SU2U1PG grp;
typedef TwoU1PG mapgrp;
#endif



int main(int argc, char ** argv)
{
    try {
        if (argc != 2) {
            if (argc != 3) {
                std::cout << "Usage      :  " << argv[0] << " <mps.h5>" << std::endl;
                std::cout << "or optional:  " << argv[0] << " <mps.h5>" << " <target Sz value>" << std::endl;
                return 1;
            }
        }

        std::string mps_in_file = argv[1];

        if (!boost::filesystem::exists(mps_in_file))
            throw std::runtime_error("input MPS " + mps_in_file + " does not exist\n");
        if (*(mps_in_file.rbegin()) == '/')
            mps_in_file.erase(mps_in_file.size()-1, 1);

        // load source MPS
        MPS<Matrix, grp> mps;
        load(mps_in_file, mps);

        // fetch parameters and modify symmetry
        storage::archive ar_in(mps_in_file + "/props.h5");
        BaseParameters parms;
        ar_in["/parameters"] >> parms;
        parms.set("init_state", "const");
#if defined(USE_SU2U1)
        parms.set("symmetry", "2u1");
#elif defined(USE_SU2U1PG)
        parms.set("symmetry", "2u1pg");
#endif

        int N = parms["nelec"];
        int TwoS = parms["spin"];

	    int target_sz = 0;
	    if(argc == 3)
            target_sz = atoi(argv[2]);

        // open file to store names of transformed MPS for preprocessing
	    std::ofstream myfile;
	    myfile.open ("mpstransform.txt");

        for (int Sz = -TwoS; Sz <= TwoS; Sz += 2)
        {
    	    // skip loop to next iteration until we hit the target Sz value 
	        if((argc == 3) && (Sz != target_sz))
	            continue;

            int Nup = (N + Sz) / 2;
            int Ndown = (N - Sz) / 2;

            parms.set("u1_total_charge1", Nup);
            parms.set("u1_total_charge2", Ndown);

            // the output MPS
            MPS<Matrix, mapgrp> mps_out = transform_mps<Matrix, grp>()(mps, Nup, Ndown);
            
            std::string mps_out_file = mps_in_file;
            std::size_t pos = mps_out_file.find(".h5");
            if (pos != mps_out_file.size())
                mps_out_file.erase(pos, 3);
            mps_out_file += "." + boost::lexical_cast<std::string>(TwoS) + "." + boost::lexical_cast<std::string>(Nup-Ndown) + ".h5";

            myfile << mps_out_file << std::endl;
	        myfile << Nup << std::endl;
	        myfile << Ndown << std::endl;
	        myfile << parms["irrep"] << std::endl;

            save(mps_out_file, mps_out);

            if (boost::filesystem::exists(mps_out_file + "/props.h5"))
                boost::filesystem::remove(mps_out_file + "/props.h5");
            boost::filesystem::copy(mps_in_file + "/props.h5", mps_out_file + "/props.h5");

            storage::archive ar_out(mps_out_file + "/props.h5", "w");
            ar_out["/parameters"] << parms;
        }

        myfile.close();
        
    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
