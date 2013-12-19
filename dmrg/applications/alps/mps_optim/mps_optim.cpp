/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include "dmrg/block_matrix/symmetry/nu1.h"
typedef NU1 grp;

#include "dmrg/block_matrix/detail/alps.hpp"
typedef alps::numeric::matrix<double>                matrix;
typedef alps::numeric::matrix<std::complex<double> > cmatrix;

#include "dmrg_sim.hpp"

#include <boost/chrono.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

int main(int argc, char ** argv)
{
    try {
        if (argc != 2 && argc != 3)
        {
            maquis::cout << "Usage: <parms>" << std::endl;
            exit(1);
        }
        std::string infile(argv[1]);
        std::string outfile(boost::replace_last_copy(infile, ".in.", ".out."));
        
        maquis::cout.precision(10);

        /// Task timer
        boost::chrono::system_clock::time_point t0, t1;
        t0 = boost::chrono::system_clock::now();

        /// Load parameters
        DmrgParameters parms;
        {
            std::ifstream param_file(infile.c_str());
            if (!param_file)
                throw std::runtime_error("Could not open parameter file.");
            alps::Parameters p; p.extract_from_xml(param_file);
            parms = DmrgParameters(p);
        }
        
        /// Match parameters of ALPS DMRG
        parms.set("ngrowsweeps", 0);
        parms.set("nsweeps", int(parms["SWEEPS"]));
        parms.set("max_bond_dimension", int(parms["MAXSTATES"]));
        parms.set("chkpfile", outfile.substr(0, outfile.find_last_of('.')) + ".chkp");
        parms.set("resultfile", outfile.substr(0, outfile.find_last_of('.')) + ".h5");
        
        /// Check which matrix to use
        if (parms["COMPLEX"]) {
            dmrg_sim<cmatrix, grp> sim(parms);
            sim.run();
        } else {
            dmrg_sim<matrix, grp> sim(parms);
            sim.run();
        }
        
        /// Final time reporting
        t1 = boost::chrono::system_clock::now();
        maquis::cout << "Task took " << boost::chrono::duration<double>(t1-t0) << " seconds." << std::endl;
    } catch (std::exception & e) {
        maquis::cerr << "Exception thrown:" << std::endl;
        maquis::cerr << e.what() << std::endl;
        exit(1);
    }
}

