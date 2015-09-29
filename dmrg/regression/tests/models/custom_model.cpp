/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2015 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/mpl/list.hpp>

#include <iterator>
#include <iostream>


using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/block_matrix/detail/alps.hpp"

#include "dmrg/utils/DmrgParameters.h"

#include "dmrg/models/custom_model.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/lattice.h"

#include "dmrg/mp_tensors/mps.h"


typedef alps::numeric::matrix<double> matrix;


struct U1System {
    
    typedef U1 grp;

    static const int L = 10;

    static Index<grp> phys_dim()
    {
        // hard-core bosons
        Index<grp> phys;
        phys.insert(std::make_pair(0, 1));
        phys.insert(std::make_pair(1, 1));
        return phys;
    }
    
    static grp::charge total_quantum_numbers()
    {
        return 5;
    }
    
    static void fill_terms_in_model(CustomModel<matrix, grp> & model_builder)
    {
        SiteOperator<matrix, grp> b, bdag;
        b.insert_block(matrix(1,1,1), 1, 0);
        bdag.insert_block(matrix(1,1,1), 0, 1);
        
        for (int i=0; i<L-1; ++i) {
            model_builder.add_bondterm(bdag, i, b,    i+1, -1.);
            model_builder.add_bondterm(b,    i, bdag, i+1, -1.);
        }
    }
};


struct TwoU1System {
    
    typedef TwoU1 grp;
    
};


//typedef boost::mpl::list<U1System, TwoU1System> test_systems;
typedef boost::mpl::list<U1System> test_systems;


BOOST_AUTO_TEST_CASE_TEMPLATE( custom_model_energy, ML, test_systems )
{
    typedef typename ML::grp grp;
    typedef typename grp::charge charge;
    
    const int L = ML::L;
    DmrgParameters parms;
    parms.set("max_bond_dimension", 40);
    parms.set("lattice_library", "coded");
    parms.set("LATTICE", "chain lattice");
    parms.set("L", L);
    parms.set("a", 1.);

    Lattice lattice(parms);

    Index<grp> phys = ML::phys_dim();
    CustomModel<matrix,grp> model_builder(phys);
    ML::fill_terms_in_model(model_builder);
    
    MPO<matrix, grp> mpo = make_mpo(lattice, model_builder.make_model());
    
    charge initc = ML::total_quantum_numbers();
    default_mps_init<matrix, grp> initializer(parms, std::vector<Index<grp> >(1, phys), initc, std::vector<int>(L,0));
    MPS<matrix, grp> mps(L, initializer);
    
    double energy = expval(mps, mpo);
    maquis::cout << "Random energy: " << energy << std::endl;
    
    BOOST_CHECK_CLOSE(energy, -5.445283602, 1e-7); // if the random state does not change, it should be the correct energy.
}

