/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN

// Unit test for integral map
#include <boost/test/included/unit_test.hpp>
#include "utils/fpcomparison.h"
#include "dmrg/block_matrix/symmetry/gsl_coupling.h"

BOOST_AUTO_TEST_CASE( Test_Wigner )
{
    // Test for calculation of Wigner 9j symbols, in gsl API (2*n)
    std::vector<std::pair<std::array<int, 9>, double> > testcases = {
        { {0,0,0, 0,0,0, 0,0,0}, 1.0 },
        { {0,0,0, 1,0,1, 1,0,1}, 0.5 },
        { {1,0,1, 0,0,0, 1,0,1}, 0.5 },
        { {0,1,1, 1,0,1, 1,1,0}, -0.25 },
        { {0,1,1, 1,0,1, 1,1,2}, 0.25 },
        { {0,1,1, 1,1,0, 1,0,1}, 0.25 },
        { {0,1,1, 1,1,0, 1,2,1}, 0.25 },
        { {0,1,1, 1,1,2, 1,0,1}, 0.25 },
        { {0,1,1, 1,2,1, 1,1,0}, 0.25 },
        { {0,1,1, 1,2,1, 1,1,2}, 0.08333333333 },
        { {0,1,1, 1,2,1, 1,2,2}, 0.0 },
        { {0,1,1, 2,0,2, 2,1,1}, -0.1666666667 },
        { {0,1,1, 2,0,2, 2,1,3}, 0.1666666667 },
        { {0,1,1, 2,1,1, 2,0,2}, 0.1666666667 },
        { {2,0,2, 2,1,1, 0,1,1}, -0.1666666667 },
        { {2,2,0, 2,1,1, 0,1,1}, 0.1666666667  },
        { {2,2,0, 2,1,1, 3,2,1}, 0 },
        { {3,1,4, 1,0,1, 3,1,3}, 0 },
        { {3,1,2, 2,1,1, 1,2,3}, -0.02777777778},
        { {3,1,2, 1,2,3, 2,1,1}, -0.02777777778},
        { {2,1,1, 3,1,2, 1,2,3}, -0.02777777778},
        { {1,2,1, 1,3,2, 2,1,3}, -0.02777777778},
        { {2,1,1, 3,2,1, 1,3,2}, -0.02777777778},
        { {7,2,7, 2,1,1, 5,1,6}, 0.03571428571 }
    };

    WignerWrapper::UseCache = true;
    WignerWrapper::fill_cache(7); // maximum value in the tests
    for (auto&& test: testcases)
    {
        auto& p = test.first;
        BOOST_CHECK_CLOSE(
            WignerWrapper::gsl_sf_coupling_9j(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8]),
            test.second,
            1e-6);
    }
}
