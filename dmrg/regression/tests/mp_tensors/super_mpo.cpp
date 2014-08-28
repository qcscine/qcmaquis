/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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

#include <cmath>
#include <iterator>
#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/block_matrix/detail/alps.hpp"

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/super_mpo.h"


typedef alps::numeric::matrix<double> matrix;


template <class SymmGroup>
void mpo1_test(Index<SymmGroup> const& phys,
               typename OPTable<matrix, SymmGroup>::tag_type op1,
               typename OPTable<matrix, SymmGroup>::tag_type op2,
               typename OPTable<matrix, SymmGroup>::tag_type op3,
               boost::shared_ptr<OPTable<matrix, SymmGroup> > op_table)
{
    maquis::cout << "TEST MPO 1" << std::endl;

    typedef typename OPTable<matrix, SymmGroup>::tag_type tag_type;
    typedef boost::tuple<std::size_t, std::size_t, tag_type, matrix::value_type> prempo_element;
    typedef std::vector<prempo_element> prempo_type;
    
    /// MPO with three sites.
    ///   - op1 - op1 - op2 - op3 -
    ///   `      `op2 - op3 - op1´
    ///    `op2 - op3 - op1´
    MPO<matrix,SymmGroup> mpo(4);
    int site = 0;
    {
        prempo_type mpo_description;
        mpo_description.push_back( prempo_element(0, 0, op1, 1.) );
        mpo_description.push_back( prempo_element(0, 1, op2, 1.) );
        
        mpo[site++] = MPOTensor<matrix, SymmGroup>(1, 2, mpo_description, op_table);
    }
    
    {
        prempo_type mpo_description;
        mpo_description.push_back( prempo_element(0, 0, op1, 1.) );
        mpo_description.push_back( prempo_element(1, 1, op3, 1.) );
        mpo_description.push_back( prempo_element(0, 2, op2, 1.) );

        mpo[site++] = MPOTensor<matrix, SymmGroup>(2, 3, mpo_description, op_table);
    }
    
    {
        prempo_type mpo_description;
        mpo_description.push_back( prempo_element(1, 1, op1, 1.) );
        mpo_description.push_back( prempo_element(2, 1, op3, 1.) );
        mpo_description.push_back( prempo_element(0, 2, op2, 1.) );

        mpo[site++] = MPOTensor<matrix, SymmGroup>(3, 3, mpo_description, op_table);
    }
    
    {
        prempo_type mpo_description;
        mpo_description.push_back( prempo_element(1, 1, op1, 1.) );
        mpo_description.push_back( prempo_element(2, 1, op3, 1.) );
        
        mpo[site++] = MPOTensor<matrix, SymmGroup>(3, 2, mpo_description, op_table);
    }
    
    
    MPS<matrix, SymmGroup> super_mps = mpo_to_smps(mpo, phys);
    std::cout << "MPS description:" << std::endl << super_mps.description();
}

template <class SymmGroup>
void mpo2_test(Index<SymmGroup> const& phys,
               typename OPTable<matrix, SymmGroup>::tag_type op1,
               typename OPTable<matrix, SymmGroup>::tag_type op2,
               typename OPTable<matrix, SymmGroup>::tag_type op3,
               boost::shared_ptr<OPTable<matrix, SymmGroup> > op_table)
{
    maquis::cout << "TEST MPO 2" << std::endl;
    
    typedef typename OPTable<matrix, SymmGroup>::tag_type tag_type;
    typedef boost::tuple<std::size_t, std::size_t, tag_type, matrix::value_type> prempo_element;
    typedef std::vector<prempo_element> prempo_type;
    
    /// MPO with three sites.
    ///   - op1 - op1 - op2 -
    ///   `      `op2 - op3 -
    ///    `op2 - op3 - op1´
    MPO<matrix,SymmGroup> mpo(3);
    int site = 0;
    {
        prempo_type mpo_description;
        mpo_description.push_back( prempo_element(0, 0, op1, 1.) );
        mpo_description.push_back( prempo_element(0, 1, op2, 1.) );
        
        mpo[site++] = MPOTensor<matrix, SymmGroup>(1, 2, mpo_description, op_table);
    }
    
    {
        prempo_type mpo_description;
        mpo_description.push_back( prempo_element(0, 0, op1, 1.) );
        mpo_description.push_back( prempo_element(1, 1, op3, 1.) );
        mpo_description.push_back( prempo_element(0, 2, op2, 1.) );

        mpo[site++] = MPOTensor<matrix, SymmGroup>(2, 3, mpo_description, op_table);
    }
    
    {
        prempo_type mpo_description;
        mpo_description.push_back( prempo_element(1, 1, op1, 1.) );
        mpo_description.push_back( prempo_element(2, 1, op3, 1.) );
        mpo_description.push_back( prempo_element(0, 2, op2, 1.) );

        mpo[site++] = MPOTensor<matrix, SymmGroup>(3, 3, mpo_description, op_table);
    }
    
    
    MPS<matrix, SymmGroup> super_mps = mpo_to_smps(mpo, phys);
    std::cout << "MPS description:" << std::endl << super_mps.description();
}

template <class SymmGroup>
void mpo3_test(Index<SymmGroup> const& phys,
               typename OPTable<matrix, SymmGroup>::tag_type op1,
               typename OPTable<matrix, SymmGroup>::tag_type op2,
               typename OPTable<matrix, SymmGroup>::tag_type op3,
               boost::shared_ptr<OPTable<matrix, SymmGroup> > op_table)
{
    maquis::cout << "TEST MPO 3" << std::endl;
    
    typedef typename OPTable<matrix, SymmGroup>::tag_type tag_type;
    typedef boost::tuple<std::size_t, std::size_t, tag_type, matrix::value_type> prempo_element;
    typedef std::vector<prempo_element> prempo_type;
    
    /// MPO with four sites.
    ///  - op1 - op3 - op1 - op1
    /// `            `           `
    ///  `            `op2 - op3 -`
    ///   `op3 - op1 ´            ´
    ///              ` op3 - op2 ´

    MPO<matrix,SymmGroup> mpo(4);
    int site = 0;
    {
        prempo_type mpo_description;
        mpo_description.push_back( prempo_element(0, 0, op1, 1.) );
        mpo_description.push_back( prempo_element(0, 1, op3, 1.) );
        
        mpo[site++] = MPOTensor<matrix, SymmGroup>(1, 2, mpo_description, op_table);
    }
    
    {
        prempo_type mpo_description;
        mpo_description.push_back( prempo_element(0, 0, op3, 1.) );
        mpo_description.push_back( prempo_element(1, 1, op1, 1.) );
        
        mpo[site++] = MPOTensor<matrix, SymmGroup>(2, 2, mpo_description, op_table);
    }
    
    {
        prempo_type mpo_description;
        // this mpo seems impossible...
        mpo_description.push_back( prempo_element(0, 0, op1, 1.) );
        mpo_description.push_back( prempo_element(0, 1, op2, 1.) );
        mpo_description.push_back( prempo_element(0, 2, op3, 1.) );
        
        mpo[site++] = MPOTensor<matrix, SymmGroup>(2, 3, mpo_description, op_table);
    }

    {
        prempo_type mpo_description;
        mpo_description.push_back( prempo_element(0, 1, op1, 1.) );
        mpo_description.push_back( prempo_element(1, 1, op3, 1.) );
        mpo_description.push_back( prempo_element(2, 1, op2, 1.) );
        
        mpo[site++] = MPOTensor<matrix, SymmGroup>(3, 3, mpo_description, op_table);
    }

    MPS<matrix, SymmGroup> super_mps = mpo_to_smps(mpo, phys);
    std::cout << "MPS description:" << std::endl << super_mps.description();
}


template <class SymmGroup>
void all_tests(Index<SymmGroup> const& phys,
               typename OPTable<matrix, SymmGroup>::tag_type op1,
               typename OPTable<matrix, SymmGroup>::tag_type op2,
               typename OPTable<matrix, SymmGroup>::tag_type op3,
               boost::shared_ptr<OPTable<matrix, SymmGroup> > op_table)
{
    
    mpo1_test(phys, op1, op2, op3, op_table);
    mpo2_test(phys, op1, op2, op3, op_table);
}

void test_none()
{
    typedef TrivialGroup SymmGroup;
    typedef OPTable<matrix, SymmGroup>::tag_type tag_type;
    maquis::cout << "TESTING NONE SYMMETRY" << std::endl;
    
    SymmGroup::charge C = SymmGroup::IdentityCharge;
    
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(C, 2));
    
    block_matrix<matrix, SymmGroup> op1, op2, op3;
    
    op1 = identity_matrix<matrix>(phys);
    std::cout << "Operator 1:" << std::endl << op1;
    
    {
        matrix tmp(2,2,0.);
        tmp(0,0) = 1.; tmp(1,1) = -1.;
        op2.insert_block(tmp, C, C);
    }
    std::cout << "Operator 2:" << std::endl << op2;
    
    {
        matrix tmp(2,2,0.); tmp(0,1) = 1.;
        op3.insert_block(tmp, C, C);
    }
    std::cout << "Operator 3:" << std::endl << op3;

    boost::shared_ptr<OPTable<matrix, SymmGroup> > op_table(new OPTable<matrix, SymmGroup>());
    tag_type op1_tag = op_table->register_op(op1);
    tag_type op2_tag = op_table->register_op(op2);
    tag_type op3_tag = op_table->register_op(op3);

    all_tests(phys, op1_tag, op2_tag, op3_tag, op_table);
}

void test_u1()
{
    typedef U1 SymmGroup;
    typedef OPTable<matrix, SymmGroup>::tag_type tag_type;
    maquis::cout << "TESTING U1 SYMMETRY" << std::endl;
    
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 2));
    
    block_matrix<matrix, SymmGroup> op1, op2, op3;
    
    op1 = identity_matrix<matrix>(phys);
    std::cout << "Operator 1:" << std::endl << op1;
    
    {
        matrix tmp(1,2,0.); tmp(0,0) = 1.;
        op2.insert_block(tmp, 0, 1);
    }
    std::cout << "Operator 2:" << std::endl << op2;
    
    {
        matrix tmp(2,1,0.); tmp(1,0) = 1.;
        op3.insert_block(tmp, 1, 0);
    }
    std::cout << "Operator 3:" << std::endl << op3;
    
    boost::shared_ptr<OPTable<matrix, SymmGroup> > op_table(new OPTable<matrix, SymmGroup>());
    tag_type op1_tag = op_table->register_op(op1);
    tag_type op2_tag = op_table->register_op(op2);
    tag_type op3_tag = op_table->register_op(op3);

    all_tests(phys, op1_tag, op2_tag, op3_tag, op_table);
}


int main() {
    try {
        test_none();
        test_u1();
    } catch (std::exception & e) {
        std::cerr << "Exception:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
