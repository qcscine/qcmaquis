/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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

#include "dmrg/sim/matrix_types.h"
#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/chem/util.h"

#include "ci_encode.hpp"

typedef TwoU1PG grp;


int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        maquis::cout << "Usage: mps2ci <mps.h5> <determinants_file> " << std::endl;
        maquis::cout << "See J. Chem. Phys. 126, 244109(2007)" << std::endl;
        exit(1);
    }

    maquis::cout.precision(10);

    std::ifstream param_file(argv[1]);
    if (!param_file) {
        maquis::cerr << "Could not open the mps." << std::endl;
        exit(1);
    }

    typedef int pos_t;

    // load the MPS
    MPS<matrix, grp> mps;
    load(argv[1], mps);

    pos_t L = mps.length();
    grp::subcharge Nup = mps[L-1].col_dim()[0].first[0];
    grp::subcharge Ndown = mps[L-1].col_dim()[0].first[1];

    BaseParameters parms;
    parms.set("site_types", chem::detail::infer_site_types(mps));

    // extract physical basis for every site from MPS
    std::vector<grp::subcharge> irreps = parms["site_types"];
    std::vector<Index<grp> > per_site, phys_dims = chem::detail::make_2u1_site_basis<matrix, grp>(L, Nup, Ndown, parms["site_types"]);
    for (pos_t q = 0; q < L; ++q)
        per_site.push_back(phys_dims[irreps[q]]);

    // load the determinants
    std::vector<std::vector<grp::charge> > determinants = parse_config<matrix, grp>(std::string(argv[2]), per_site);
    // printout the determinants
    for (pos_t q = 0;q < determinants.size(); ++q){
       for (pos_t p = 0; p < L; ++p){
           std::cout << determinants[q][p];
       }
       std::cout << std::endl;
    }

    // compute the CI coefficients for all determinants in the input
    int i = 1;
    for (std::vector< std::vector<grp::charge> >::iterator it = determinants.begin();
        it != determinants.end(); ++it){
        maquis::cout << "CI coefficient of det " << i <<": " << extract_coefficient(mps, *it) << std::endl;
        i++;
        }

    maquis::cout << std::endl;

    return 0;
}

