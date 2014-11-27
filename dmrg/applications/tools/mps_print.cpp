/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

//#include <cmath>
#include <iostream>
#include <fstream>
//#include <iomanip>
//#include <sys/time.h>
//#include <sys/stat.h>

//#include <boost/lambda/lambda.hpp>

//using std::cerr;

//#include <alps/numeric/matrix.hpp>
//#include <alps/numeric/matrix/algorithms.hpp>
#include "dmrg/block_matrix/detail/alps.hpp"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpstensor.h"
//#include "dmrg/mp_tensors/mpo.h"
//#include "dmrg/mp_tensors/boundary.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
//#include "dmrg/mp_tensors/contractions.h"

typedef alps::numeric::matrix<double> matrix;

typedef TwoU1LPG grp;

int main(int argc, char ** argv)
{
    try {
        if (argc < 2) {
            std::cout << "Usage: " << argv[0] << " <mps.h5>" << std::endl;
            return 1;
        }
        //cout.precision(5);

        MPS<matrix, grp> mps;
        load(argv[1], mps);
        size_t L = mps.size();

        std::vector<int> site_irreps;
        for (int i=0; i < L; ++i)
            site_irreps.push_back(mps[i].site_dim()[0].first[2]);

        std::cout << "site irreps: ";
        std::copy(site_irreps.begin(), site_irreps.end(), std::ostream_iterator<int>(std::cout, " "));
        std::cout << std::endl;

        std::cout << "Printing left paired mps:" << std::endl;
        for (int i=0; i < L; ++i) {
            std::cout << "################## SITE " << i << " ##################\n" << std::endl;
            MPSTensor<matrix, grp> tmp(mps[i]);
            tmp.make_left_paired();
            std::cout << tmp;
        }

        bool print_flag = false;
        if(print_flag) {
            std::ofstream myfile;
            myfile.open ("mps_leftpaired");
            myfile << "***** Left paired MPS *****\n" << "Number of sites: " << L << std::endl;
            myfile << "Site irreps are: ";
            for(int i=0; i < L; ++i) {
                myfile << site_irreps[i] << " ";
            }
            myfile << std::endl << std::endl;
            for (int i=0; i < L; ++i) {
                myfile << "################## SITE " << i << " ##################\n\n";
                MPSTensor<matrix, grp> tmp(mps[i]);
                tmp.make_left_paired();
                myfile << tmp;
            }
            myfile.close();
        }

        //Check <psi|psi> == 1
        //double norm_expectation = overlap(mps,mps);
        double norm_expectation = norm(mps);

        std::cout << "Norm of the ground state wave function:\n<psi|psi> = " << norm_expectation << std::endl;


    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
