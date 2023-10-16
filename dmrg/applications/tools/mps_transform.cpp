/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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

#include "dmrg/sim/matrix_types.h"
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
        MPS<matrix, grp> mps;
        load(mps_in_file, mps);

        // fetch parameters and modify symmetry
        storage::archive ar_in(mps_in_file + "/props.h5");
        BaseParameters parms;
        ar_in["/parameters"] >> parms;
        parms.set("init_type", "const");
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
            MPS<matrix, mapgrp> mps_out = transform_mps<matrix, grp>()(mps, Nup, Ndown);
            
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
