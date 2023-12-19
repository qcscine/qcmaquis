/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "MPSOverlapClass.h"
#include "dmrg/sim/matrix_types.h"

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

/** @brief Utility needed to calculate the overlap between two MPSs */

int main(int argc, char ** argv)
{
    try {
        if (argc != 3) {
            std::cout << "Usage: " << argv[0] << " <mps1.h5> <mps2.h5>" << std::endl;
            return 1;
        }
        // Creates the overlap calculator object and calculates the overlap
        MPSOverlapClass<matrix, grp> overlapCalculator(argv[1], argv[2]);
        overlapCalculator.printOverlap();

        // OLD CODE - TO CHECK IF IT'S ACTUALLY USEFUL
        //operator_selector<matrix, grp>::type ident;
        //for (int i=0; i<mps1.site_dim(0).size(); ++i)
        //    ident.insert_block(matrix::identity_matrix(mps1.site_dim(0)[i].second),
        //                       mps1.site_dim(0)[i].first, mps1.site_dim(0)[i].first);
        //
        //MPO<matrix, grp> mpo;
        //
        //MPOTensor<matrix, grp> mpot;
        //mpot.set(0,0, ident);
        //mpo = MPO<matrix, grp>(mps1.length());
        //for (int p=0; p<mps1.length(); ++p)
        //    mpo[p] = mpot;
        //
        //std::cout << "<mps1 | 1 | mps2> = " << expval(mps1, mps2, mpo) << std::endl;
        //std::cout << "<mps2 | 1 | mps1> = " << expval(mps2, mps1, mpo) << std::endl;
    }
    catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
