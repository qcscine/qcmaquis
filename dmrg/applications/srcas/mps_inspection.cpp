/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifdef USE_AMBIENT
    #include "dmrg/block_matrix/detail/ambient.hpp"
    typedef ambient::numeric::tiles<ambient::numeric::matrix<double> > Matrix;
#else
    #include "dmrg/block_matrix/detail/alps.hpp"
    typedef alps::numeric::matrix<double> Matrix;
#endif

#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/boundary.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mps.h"

#include "dmrg/mp_tensors/mps_initializers.h"

#include "ci_encode.hpp"
#include "ci_deas.hpp"

typedef TwoU1 grp; 


template<class Matrix>
std::pair<std::vector<typename Matrix::value_type>, std::vector<std::vector<typename TwoU1::charge> > >
parse_config(std::string file)
{
    std::ifstream config_file;
    config_file.open(file.c_str());

    TwoU1::charge A(1), B(0), C(0), D(0);
    B[0]=1; C[1]=1;

    std::vector<typename Matrix::value_type> coefficients;
    std::vector<std::vector<typename TwoU1::charge> > configs;

    for (std::string line; std::getline(config_file, line); ) {
        //boost::char_separator<char> sep(" ");
        std::vector<std::string> det_coeff;
        boost::split(det_coeff, line, boost::is_any_of(" "));
        
        std::string det = det_coeff[0];
        typename Matrix::value_type coeff = boost::lexical_cast<typename Matrix::value_type>(det_coeff[1]); 
                
        coefficients.push_back(coeff);

        std::vector<typename TwoU1::charge> tmp;
        for (std::size_t i = 0; i < det.size(); ++i) {
            int occ = boost::lexical_cast<int>(det[i]); 
            switch(occ) {
                case 4:
                    tmp.push_back(A);
                    break;
                case 3:
                    tmp.push_back(B);
                    break;
                case 2:
                    tmp.push_back(C);
                    break;
                case 1:
                    tmp.push_back(D);
                    break;
            }
        }

        configs.push_back(tmp);
        
    }
    
    return std::make_pair(coefficients, configs);
}

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        maquis::cout << "Usage: <parms>" << std::endl;
        exit(1);
    }
    
    maquis::cout.precision(10);

    std::ifstream param_file(argv[1]);
    if (!param_file) {
        maquis::cerr << "Could not open parameter file." << std::endl;
        exit(1);
    }
    DmrgParameters parms(param_file);
    
    typedef int pos_t;

    // *** Data setup *** //
    Index<grp> phys;
    grp::charge A(1), B(0), C(0), D(0);
    B[0]=1; C[1]=1;
    phys.insert(std::make_pair(A, 1));
    phys.insert(std::make_pair(B, 1));
    phys.insert(std::make_pair(C, 1));
    phys.insert(std::make_pair(D, 1));

    grp::charge right_end;
    right_end[0] = parms["u1_total_charge1"];
    right_end[1] = parms["u1_total_charge2"];
    pos_t L = pos_t(parms["L"]);
    maquis::cout << "Total charge: " << right_end << std::endl;
    
    std::vector<Index<grp> > phys_dims(L, phys);
    std::vector<int> site_types(L, 0);
    
    default_mps_init<Matrix, grp> di(parms, phys_dims, right_end, site_types);
    MPS<Matrix, grp> mps1(L, di), mps2(L, di);
    for (pos_t p = 0; p < L; ++p)
        mps1[p].multiply_by_scalar(0.);

    // *** Construct MPS from CI coefficients *** //
    /*
    std::vector< std::vector<grp::charge> > configs;
    std::vector< typename Matrix::value_type> coeffs;

    boost::tie(coeffs, configs) = parse_config<Matrix>("config");

    for (typename std::vector< std::vector<grp::charge> >::iterator it = configs.begin();
            it != configs.end(); ++it)
        set_coefficient(mps1, *it, coeffs[std::distance(configs.begin(), it)]);

    for (typename std::vector< std::vector<grp::charge> >::iterator it = configs.begin();
            it != configs.end(); ++it)
        maquis::cout << "CI coefficient: " << extract_coefficient(mps1, *it) << std::endl;

    maquis::cout << std::endl;
    

    for (pos_t p = 0; p < L; ++p)
        maquis::cout << mps1[p] << std::endl; 
    */

    std::vector<std::string> environment = parse_config_strings("config");

    std::vector< std::map<grp::charge, std::set<std::string> > > data;
    data = arrange_configs(environment);

    display_environment(data);
        
    return 0;
}
