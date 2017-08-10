/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2014-2014 by Yingjin Ma       <yingjin.ma@phys.ethz.ch>
 *               2017-2017 by Alberto Baiardi  <alberto.baiardi@sns.it>
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

#include <cmath>
#include <iterator>
#include <iostream>
#include <string>
#include <sys/time.h>
#include <sys/stat.h>

#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/block_matrix/detail/alps.hpp"
typedef alps::numeric::matrix<double> matrix;

#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"

#include "ci_encode.hpp"
#include "sampling.hpp"

typedef TrivialGroup grp;

//
// Small function to read determinants from a input file
// -----------------------------------------------------

template<class Matrix, class SymmGroup>
std::vector<std::vector<typename SymmGroup::charge> >
parse_config(std::string file, std::vector<Index<SymmGroup> > const & site_dims, std::vector<int> order)
{
    // Initialization
    std::ifstream config_file;
    config_file.open(file.c_str());
    typename SymmGroup::charge A(1), B(0), C(0), D(0);
    B[0]=1;
    C[1]=1;
    std::vector<std::vector<typename SymmGroup::charge> > configs;
    // Parsing of the determinant
    for (std::string line; std::getline(config_file, line); ) {
        std::vector<std::string> det_coeff;
        boost::split(det_coeff, line, boost::is_any_of(" "));
        std::string det = det_coeff[0];
        if (det.size() != site_dims.size())
            throw std::runtime_error("The determinant length doesn't match the mps length\n");
        std::vector<typename SymmGroup::charge> tmp;
        for (std::size_t i = 0; i < det.size(); ++i) {
            int occ = boost::lexical_cast<int>(det[order[i]-1]); 
            switch(occ) {
                case 4:
                    tmp.push_back(site_dims[i][0].first); // doubly occ
                    break;
                case 3:
                    tmp.push_back(site_dims[i][1].first); // up
                    break;
                case 2:
                    tmp.push_back(site_dims[i][2].first); // down 
                    break;
                case 1:
                    tmp.push_back(site_dims[i][3].first); // empty
                    break;
            }
        }
        configs.push_back(tmp);
    }
    return configs;
}

// +------------+
//  MAIN PROGRAM
// +------------+
//
// Program that takes in input a file contaning the info for an MPS and reconstruct
// the coefficients of a CI expansion 

int main(int argc, char ** argv)
{
    // If less than three arguments have been provided, just prints info
    if (argc < 3) {
        maquis::cout << "srcas <mps.h5> <determinants_file> <CI_threshold> <CI_completeness> ( <Nsamples>  <determinants_reservoir> )" << std::endl;
        maquis::cout << "See J. Chem. Phys. 134, 224101 (2011) for details" << std::endl;
        exit(1);
    }
    maquis::cout.precision(10);
    std::ifstream param_file(argv[1]);
    if (!param_file) {
        maquis::cerr << "Could not open the mps." << std::endl;
        exit(1);
    }
    typedef int pos_t;
    //
    // LOADING OF INPUT DATA
    // ---------------------
    // Load the MPS
    MPS<matrix, grp> mps;
    load(argv[1], mps);
    pos_t L = mps.length();
    // Loading the archive:
    storage::archive ar_in(argv[1]+std::string("/props.h5"));
    // Loading the parameters
    DmrgParameters parms;
    ar_in["/parameters"] >> parms;
    // Reads ordering from parms:
    std::vector<pos_t> order(mps.length());
    order = parms["orbital_order"].as<std::vector<pos_t> >();
    maquis::cout << " orbital_order " << std::endl;
    for (pos_t p = 0; p < L; ++p)
        maquis::cout << " " << order[p] << " "; // Here should be mps[order[p]]
    std::cout << std::endl; 
    // Extract physical basis for every site from MPS
    std::vector<Index<grp> > phys_dims;
    for (pos_t p = 0; p < L; ++p)
        phys_dims.push_back(mps[p].site_dim());
    // Load the determinants
    std::vector<std::vector<typename grp::charge> > determinants = parse_config<matrix, grp>(std::string(argv[2]), phys_dims, order);
    // printout the determinants
    for (pos_t p = 0; p < L; ++p)
        std::cout << determinants[0][p];
    std::cout << std::endl;
    if(argc == 3){
        // compute the CI coefficients for all determinants in the input
        for (typename std::vector< std::vector<grp::charge> >::iterator it = determinants.begin(); it != determinants.end(); ++it)
            maquis::cout << "CI coefficient: " << extract_coefficient(mps, *it) << std::endl;
    } else {
        // Sets default input parameters
        int     Nsamples      = 1000 ;
        double  CI_threshold  = 0.01 ;
        double  COM_threshold = 0.01 ;
        // Changes if parameters have been provided in input
        if(argc > 3){
            CI_threshold  = atof(argv[3]);
        } else if(argc > 4) {
            COM_threshold = atof(argv[4]);
        } else if(argc > 5) {
             Nsamples = atoi(argv[5]);
        }
        // determinants reservoir | optional
        std::vector<std::vector<typename grp::charge> > determinants_mclr;
        if(argc == 7){
            determinants_mclr = parse_config<matrix, grp>(std::string(argv[6]), phys_dims, order);
        } else {
            determinants_mclr = determinants;
        }
        // This is used for determinants reservoir
        typedef std::map<std::vector<typename grp::charge>, double> Hash_Map_with_value ;
        Hash_Map_with_value hash_value ;
        // determinants index in determinants reservoir
        typedef std::map<std::vector<typename grp::charge>, long>   Hash_Map_with_index ;
        Hash_Map_with_index hash_index;
        // SR-CAS -- Sampling Reconstructed CAS determinants
        Sampling().generate_dets <std::vector<std::vector<typename grp::charge> >, MPS<matrix, grp>, Hash_Map_with_value, Hash_Map_with_index, std::vector<Index<grp> >, std::vector<typename grp::charge>, grp >
                (determinants, determinants_mclr, mps, hash_value, hash_index, phys_dims, L, Nsamples, CI_threshold, COM_threshold);
    }
    maquis::cout << std::endl;
    return 0;
}

