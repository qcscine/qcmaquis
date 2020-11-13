/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2018         Leon Freitag <lefreita@ethz.ch>
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

#define BOOST_TEST_MAIN

// Unit test for integral map
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
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
