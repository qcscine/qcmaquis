/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
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

#include "dmrg/utils/DmrgParameters2.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_initializers.h"

#include "ci_encode.hpp"

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
    if (argc != 3)
    {
        maquis::cout << "Usage: <parms> <model_parms>" << std::endl;
        exit(1);
    }
    
    maquis::cout.precision(10);

    std::ifstream param_file(argv[1]);
    if (!param_file) {
        maquis::cerr << "Could not open parameter file." << std::endl;
        exit(1);
    }
    DmrgParameters parms(param_file);
    
    std::ifstream model_file(argv[2]);
    if (!model_file) {
        maquis::cerr << "Could not open model file." << std::endl;
        exit(1);
    }
    ModelParameters model(model_file);
    
    typedef unsigned pos_t;

    // *** Data setup *** //
    Index<grp> phys;
    grp::charge A(1), B(0), C(0), D(0);
    B[0]=1; C[1]=1;
    phys.insert(std::make_pair(A, 1));
    phys.insert(std::make_pair(B, 1));
    phys.insert(std::make_pair(C, 1));
    phys.insert(std::make_pair(D, 1));

    grp::charge right_end;
    right_end[0] = model["u1_total_charge1"];
    right_end[1] = model["u1_total_charge2"];
    pos_t L = pos_t(model["L"]);
    maquis::cout << "Total charge: " << right_end << std::endl;
    
    MPS<Matrix, grp> mps1(L), mps2(L);
    default_mps_init<Matrix, grp> di;
    di(mps1, 256, phys, right_end); 
    for (pos_t p = 0; p < L; ++p)
        mps1[p].multiply_by_scalar(0.);

    // ****************** //
    //grp::charge tmp[] = {A,A,D,D,D,D,A,A,D};
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
        
    return 0;
}
