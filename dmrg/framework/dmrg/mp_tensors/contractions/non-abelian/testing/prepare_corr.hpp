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

#ifndef CONTRACTIONS_SU2_PREPARE_OPS_HPP
#define CONTRACTIONS_SU2_PREPARE_OPS_HPP

#include <vector>

#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/block_matrix/block_matrix.h"

#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/mpo.h"

namespace SU2 {

    template<class SymmGroup>
    Index<SymmGroup> SU2phys(int Ilocal)
    {
        typename SymmGroup::charge A(0), B(0), C(0), D(0);
        A[0] = 2; // 200
        B[0] = 1; B[1] =  1; B[2] = Ilocal; // 11I
        C[0] = 1; C[1] = -1; C[2] = Ilocal; // 1-1I
        // D = 000

        Index<SymmGroup> phys_i;
        phys_i.insert(std::make_pair(A, 1));
        phys_i.insert(std::make_pair(B, 1));
        phys_i.insert(std::make_pair(C, 1));
        phys_i.insert(std::make_pair(D, 1));

        return phys_i;
    }

    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup> make_1rdm_term(int i, int j, std::vector<int> site_irreps)
    {
        //boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
        MPO<Matrix, SymmGroup> ret(site_irreps.size());
        for (int p=0; p<site_irreps.size(); ++p)
        {
            typedef tag_detail::tag_type tag_type;
            typename SymmGroup::charge A(0), B(0), C(0), D(0);
            A[0] = 2; // 200
            B[0] = 1; B[1] =  1; B[2] = site_irreps[p]; // 11I
            C[0] = 1; C[1] = -1; C[2] = site_irreps[p]; // 1-1I
            // D = 000

            block_matrix<Matrix, SymmGroup> identity;
            identity.insert_block(Matrix(1,1,1), A, A);
            identity.insert_block(Matrix(1,1,1), B, B);
            identity.insert_block(Matrix(1,1,1), C, C);
            identity.insert_block(Matrix(1,1,1), D, D);

            block_matrix<Matrix, SymmGroup> fill;
            fill.insert_block(Matrix(1,1,1), A, A);
            fill.insert_block(Matrix(1,1,1), D, D);  // c^dag * c
            fill.insert_block(Matrix(1,1,-1), B, B); // -1
            fill.insert_block(Matrix(1,1,-1), C, C); // -1
            fill.insert_block(Matrix(1,1,-1), B, C);  // -1
            fill.insert_block(Matrix(1,1,-1), C, B);  // -1

            block_matrix<Matrix, SymmGroup> destroy;
            destroy.twoS = 1; destroy.twoSaction = 1;
            destroy.insert_block(Matrix(1,1,1), A, B);        // 1 
            destroy.insert_block(Matrix(1,1,1), A, C);        // 1
            destroy.insert_block(Matrix(1,1,-sqrt(2.)), B, D); // 1
            destroy.insert_block(Matrix(1,1,-sqrt(2.)), C, D); // 1

            block_matrix<Matrix, SymmGroup> create;
            create.twoS = 1; create.twoSaction = -1;
            create.insert_block(Matrix(1,1,2.), B, A);       // 1
            create.insert_block(Matrix(1,1,2.), C, A);       // 1
            create.insert_block(Matrix(1,1,-sqrt(2.)), D, B); // 1
            create.insert_block(Matrix(1,1,-sqrt(2.)), D, C); // 1

            //tag_type ident = tag_handler->register_op(identity, tag_detail::bosonic);
            MPOTensor<Matrix, SymmGroup> op(1,1);

            if (p == i) {
                op.set(0,0,destroy,-1.0);
            }
            else if (p == j)
                op.set(0,0,create);
            else if ( i < p && p < j)
                op.set(0,0,fill, 1.0);
            else 
                op.set(0,0,identity, 1.0);

            ret[p] = op;
        }
        return ret;
    }

    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup> make_1rdm_termB(int i, int j, std::vector<int> site_irreps)
    {   // represents a cdag_i c_j operator term
        //boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
        MPO<Matrix, SymmGroup> ret(site_irreps.size());
        for (int p=0; p<site_irreps.size(); ++p)
        {
            typedef tag_detail::tag_type tag_type;
            typename SymmGroup::charge A(0), B(0), C(0), D(0);
            A[0] = 2; // 200
            B[0] = 1; B[1] =  1; B[2] = site_irreps[p]; // 11I
            C[0] = 1; C[1] = -1; C[2] = site_irreps[p]; // 1-1I
            // D = 000

            block_matrix<Matrix, SymmGroup> identity;
            identity.insert_block(Matrix(1,1,1), A, A);
            identity.insert_block(Matrix(1,1,1), B, B);
            identity.insert_block(Matrix(1,1,1), C, C);
            identity.insert_block(Matrix(1,1,1), D, D);

            block_matrix<Matrix, SymmGroup> fill;
            fill.insert_block(Matrix(1,1,1), A, A);
            fill.insert_block(Matrix(1,1,1), D, D); 
            fill.insert_block(Matrix(1,1,-1), B, B);
            fill.insert_block(Matrix(1,1,-1), C, C);
            fill.insert_block(Matrix(1,1,-1), B, C);
            fill.insert_block(Matrix(1,1,-1), C, B);

            block_matrix<Matrix, SymmGroup> create;
            create.twoS = 1; create.twoSaction = 1;
            create.insert_block(Matrix(1,1,2.), B, A);      
            create.insert_block(Matrix(1,1,2.), C, A);      
            create.insert_block(Matrix(1,1,sqrt(2.)), D, B);
            create.insert_block(Matrix(1,1,sqrt(2.)), D, C);

            block_matrix<Matrix, SymmGroup> destroy;
            destroy.twoS = 1; destroy.twoSaction = -1;
            destroy.insert_block(Matrix(1,1,1), A, B);        
            destroy.insert_block(Matrix(1,1,1), A, C);       
            destroy.insert_block(Matrix(1,1,sqrt(2.)), B, D);
            destroy.insert_block(Matrix(1,1,sqrt(2.)), C, D);

            //tag_type ident = tag_handler->register_op(identity, tag_detail::bosonic);
            MPOTensor<Matrix, SymmGroup> op(1,1);

            if (p == i)
                op.set(0,0,create,1.0);
            else if (p == j)
                op.set(0,0,destroy);
            else if ( i < p && p < j)
                op.set(0,0,fill, 1.0);
            else 
                op.set(0,0,identity, 1.0);

            ret[p] = op;
        }
        return ret;
    }

    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup> make_custom(int i, int j, std::vector<int> site_irreps)
    {
        typedef tag_detail::tag_type tag_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
        typedef typename Matrix::value_type value_type;
        typedef boost::tuple<std::size_t, std::size_t, tag_type, value_type> tag_block;
        boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler(new TagHandler<Matrix, SymmGroup>());

        MPO<Matrix, SymmGroup> ret(site_irreps.size());
        for (int p=0; p<site_irreps.size(); ++p)
        {
            typedef tag_detail::tag_type tag_type;
            typename SymmGroup::charge A(0), B(0), C(0), D(0);
            A[0] = 2; // 200
            B[0] = 1; B[1] =  1; B[2] = site_irreps[p]; // 11I
            C[0] = 1; C[1] = -1; C[2] = site_irreps[p]; // 1-1I
            // D = 000

            block_matrix<Matrix, SymmGroup> identity_op;
            identity_op.insert_block(Matrix(1,1,1), A, A);
            identity_op.insert_block(Matrix(1,1,1), B, B);
            identity_op.insert_block(Matrix(1,1,1), C, C);
            identity_op.insert_block(Matrix(1,1,1), D, D);

            block_matrix<Matrix, SymmGroup> fill_op;
            fill_op.insert_block(Matrix(1,1,1), A, A);
            fill_op.insert_block(Matrix(1,1,1), D, D);  // c^dag * c
            fill_op.insert_block(Matrix(1,1,-1), B, B); // -1
            fill_op.insert_block(Matrix(1,1,-1), C, C); // -1
            fill_op.insert_block(Matrix(1,1,1), B, C);  // -1
            fill_op.insert_block(Matrix(1,1,1), C, B);  // -1

            block_matrix<Matrix, SymmGroup> create_op;
            create_op.insert_block(Matrix(1,1,-2.), B, A);       // 1
            create_op.insert_block(Matrix(1,1,2.), C, A);       // 1
            create_op.insert_block(Matrix(1,1,-sqrt(2.)), D, B); // 1
            create_op.insert_block(Matrix(1,1,sqrt(2.)), D, C); // 1

            block_matrix<Matrix, SymmGroup> destroy_op;
            destroy_op.insert_block(Matrix(1,1,1), A, B);        // 1 
            destroy_op.insert_block(Matrix(1,1,-1), A, C);        // 1
            destroy_op.insert_block(Matrix(1,1,sqrt(2.)), B, D); // 1
            destroy_op.insert_block(Matrix(1,1,-sqrt(2.)), C, D); // 1

            tag_type identity = tag_handler->register_op(identity_op, tag_detail::bosonic);
            tag_type fill = tag_handler->register_op(fill_op, tag_detail::bosonic);
            tag_type create = tag_handler->register_op(create_op, tag_detail::bosonic);
            tag_type destroy = tag_handler->register_op(destroy_op, tag_detail::bosonic);

            std::vector<tag_block> pre_tensor;

            if (p == i) {
                pre_tensor.push_back(tag_block(0,0,destroy, 1.0));
            }
            else if (p == j) {
                pre_tensor.push_back(tag_block(0,0,create, 1.0));
                pre_tensor.push_back(tag_block(0,1,fill, 1.0));
            }
            else if (p == j+1) {
                pre_tensor.push_back(tag_block(0,0,identity, 1.0));
                pre_tensor.push_back(tag_block(1,0,create, 1.0));
            }
            else if ( i < p && p < j){
                pre_tensor.push_back(tag_block(0,0,fill, 1.0));
            }
            else 
                pre_tensor.push_back(tag_block(0,0,identity, 1.0));

            std::pair<index_type, index_type> rcd = ::generate_mpo::rcdim(pre_tensor);
            ret[p] = MPOTensor<Matrix, SymmGroup>(rcd.first, rcd.second, pre_tensor, tag_handler->get_operator_table());
        }
        return ret;
    }

    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup> make_2rdm_term(int i, int j, int k, int l, std::vector<int> site_irreps)
    {
        MPO<Matrix, SymmGroup> ret(site_irreps.size());
        for (int p=0; p<site_irreps.size(); ++p)
        {
            typedef tag_detail::tag_type tag_type;
            typename SymmGroup::charge A(0), B(0), C(0), D(0);
            A[0] = 2; // 200
            B[0] = 1; B[1] =  1; B[2] = site_irreps[p]; // 11I
            C[0] = 1; C[1] = -1; C[2] = site_irreps[p]; // 1-1I
            // D = 000

            block_matrix<Matrix, SymmGroup> identity;
            identity.insert_block(Matrix(1,1,1), A, A);
            identity.insert_block(Matrix(1,1,1), B, B);
            identity.insert_block(Matrix(1,1,1), C, C);
            identity.insert_block(Matrix(1,1,1), D, D);

            block_matrix<Matrix, SymmGroup> fill;
            fill.insert_block(Matrix(1,1,1), A, A);
            fill.insert_block(Matrix(1,1,1), D, D);  // c^dag * c
            fill.insert_block(Matrix(1,1,-1), B, B); // -1
            fill.insert_block(Matrix(1,1,-1), C, C); // -1
            fill.insert_block(Matrix(1,1,1), B, C);  // -1
            fill.insert_block(Matrix(1,1,1), C, B);  // -1

            block_matrix<Matrix, SymmGroup> createS1;
            createS1.twoS = 1; createS1.twoSaction = 1;
            createS1.insert_block(Matrix(1,1,std::sqrt(2.)), B, A);      
            createS1.insert_block(Matrix(1,1,std::sqrt(2.)), C, A);      
            createS1.insert_block(Matrix(1,1,1), D, B);
            createS1.insert_block(Matrix(1,1,1), D, C);

            block_matrix<Matrix, SymmGroup> createS0;
            createS0.twoS = 1; createS0.twoSaction = -1;
            createS0.insert_block(Matrix(1,1,2.), B, A);      
            createS0.insert_block(Matrix(1,1,2.), C, A);      
            createS0.insert_block(Matrix(1,1,-std::sqrt(2.)), D, B);
            createS0.insert_block(Matrix(1,1,-std::sqrt(2.)), D, C);

            block_matrix<Matrix, SymmGroup> createS2;
            createS2.twoS = 1; createS2.twoSaction = 1;
            createS2.insert_block(Matrix(1,1,-std::sqrt(2.)), B, A);      
            createS2.insert_block(Matrix(1,1,-std::sqrt(2.)), C, A);      
            createS2.insert_block(Matrix(1,1,1), D, B);
            createS2.insert_block(Matrix(1,1,1), D, C);

            block_matrix<Matrix, SymmGroup> destroyS0;
            destroyS0.twoS = 1; destroyS0.twoSaction = 1;
            destroyS0.insert_block(Matrix(1,1,1), A, B);      
            destroyS0.insert_block(Matrix(1,1,1), A, C);     
            destroyS0.insert_block(Matrix(1,1,-std::sqrt(2.)), B, D); 
            destroyS0.insert_block(Matrix(1,1,-std::sqrt(2.)), C, D);

            double root32 = std::sqrt(3./2.);
            block_matrix<Matrix, SymmGroup> destroyS2;
            destroyS2.twoS = 1; destroyS2.twoSaction = -1;
            destroyS2.insert_block(Matrix(1,1,root32), A, B);        
            destroyS2.insert_block(Matrix(1,1,root32), A, C);       
            destroyS2.insert_block(Matrix(1,1,-root32 * sqrt(2.)), B, D); 
            destroyS2.insert_block(Matrix(1,1,-root32 * sqrt(2.)), C, D);

            block_matrix<Matrix, SymmGroup> destroyS1;
            destroyS1.twoS = 1; destroyS1.twoSaction = -1;
            destroyS1.insert_block(Matrix(1,1,std::sqrt(2.)), A, B);      
            destroyS1.insert_block(Matrix(1,1,std::sqrt(2.)), A, C);     
            destroyS1.insert_block(Matrix(1,1,2.), B, D); 
            destroyS1.insert_block(Matrix(1,1,2.), C, D);

            //tag_type ident = tag_handler->register_op(identity, tag_detail::bosonic);
            MPOTensor<Matrix, SymmGroup> op(1,1);

            if (p == i) {
                op.set(0,0,createS1, 1.0);
            }
            else if (p == j) {
                op = MPOTensor<Matrix, SymmGroup>(1,2);
                op.set(0,0,createS0, 0.5);
                op.set(0,1,createS2, -1.0);
            }
            else if (p == k) {
                op = MPOTensor<Matrix, SymmGroup>(2,1);
                op.set(0,0,destroyS0, 1.0);
                op.set(1,0,destroyS2, 1.0);
            }
            else if (p == l) {
                op.set(0,0,destroyS1, 1.0);
            }
            //else if ( i < p && p < j)
            //    op.set(0,0,fill, 1.0);
            else  {
                op.set(0,0,identity, 1.0);
            }

            ret[p] = op;
        }
        return ret;
    }

    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup> make_count(int i, std::vector<int> site_irreps)
    {
        //boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
        MPO<Matrix, SymmGroup> ret(site_irreps.size());
        for (int p=0; p<site_irreps.size(); ++p)
        {
            typedef tag_detail::tag_type tag_type;
            typename SymmGroup::charge A(0), B(0), C(0), D(0);
            A[0] = 2; // 200
            B[0] = 1; B[1] =  1; B[2] = site_irreps[p]; // 11I
            C[0] = 1; C[1] = -1; C[2] = site_irreps[p]; // 1-1I
            // D = 000

            block_matrix<Matrix, SymmGroup> ident;
            ident.insert_block(Matrix(1,1,1), A, A);
            ident.insert_block(Matrix(1,1,1), B, B);
            ident.insert_block(Matrix(1,1,1), C, C);
            ident.insert_block(Matrix(1,1,1), D, D);

            block_matrix<Matrix, SymmGroup> count;
            count.insert_block(Matrix(1,1,2), A, A);
            count.insert_block(Matrix(1,1,1), B, B);
            count.insert_block(Matrix(1,1,1), C, C);

            MPOTensor<Matrix, SymmGroup> op(1,1);
            if (p==i)
                op.set(0,0,count, 1.0);
            else
                op.set(0,0,ident, 1.0);


            ret[p] = op;
        }
        return ret;
    }

    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup> make_e2d_d2e(int i, int j, std::vector<int> site_irreps)
    {
        typedef block_matrix<Matrix, SymmGroup> op_t;
        //boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
        MPO<Matrix, SymmGroup> ret(site_irreps.size());
        for (int p=0; p<site_irreps.size(); ++p)
        {
            typedef tag_detail::tag_type tag_type;
            typename SymmGroup::charge A(0), B(0), C(0), D(0);
            A[0] = 2; // 200
            B[0] = 1; B[1] =  1; B[2] = site_irreps[p]; // 11I
            C[0] = 1; C[1] = -1; C[2] = site_irreps[p]; // 1-1I
            // D = 000

            block_matrix<Matrix, SymmGroup> ident;
            ident.insert_block(Matrix(1,1,1), A, A);
            ident.insert_block(Matrix(1,1,1), B, B);
            ident.insert_block(Matrix(1,1,1), C, C);
            ident.insert_block(Matrix(1,1,1), D, D);

            op_t e2d_op;
            e2d_op.insert_block(Matrix(1,1,1), D, A);

            op_t d2e_op;
            d2e_op.insert_block(Matrix(1,1,1), A, D);

            MPOTensor<Matrix, SymmGroup> op(1,1);
            if (p==i)
                op.set(0,0, e2d_op, 1.0);
            else if (p==j)
                op.set(0,0, d2e_op, 1.0);
            else
                op.set(0,0, ident, 1.0);

            ret[p] = op;
        }
        return ret;
    }

    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup> make_flip(int i, int j, std::vector<double> & values, std::vector<int> site_irreps)
    {
        typedef block_matrix<Matrix, SymmGroup> op_t;
        MPO<Matrix, SymmGroup> ret(site_irreps.size());

        for (int p=0; p<site_irreps.size(); ++p)
        {
            typedef tag_detail::tag_type tag_type;
            typename SymmGroup::charge A(0), B(0), C(0), D(0);
            A[0] = 2; // 200
            B[0] = 1; B[1] =  1; B[2] = site_irreps[p]; // 11I
            C[0] = 1; C[1] = -1; C[2] = site_irreps[p]; // 1-1I
            // D = 000

            block_matrix<Matrix, SymmGroup> ident;
            ident.twoS = 0; ident.twoSaction = 0;
            ident.insert_block(Matrix(1,1,1), A, A);
            ident.insert_block(Matrix(1,1,1), B, B);
            ident.insert_block(Matrix(1,1,1), C, C);
            ident.insert_block(Matrix(1,1,1), D, D);

            block_matrix<Matrix, SymmGroup> fill;
            fill.twoS = 0; fill.twoSaction = 0;
            fill.insert_block(Matrix(1,1,1), A, A);
            fill.insert_block(Matrix(1,1,1), B, B);
            fill.insert_block(Matrix(1,1,1), C, C);
            fill.insert_block(Matrix(1,1,1), B, C);
            fill.insert_block(Matrix(1,1,1), C, B);
            fill.insert_block(Matrix(1,1,1), D, D);

            op_t flip1; flip1.twoS = 2; flip1.twoSaction = 2;
            flip1.insert_block(Matrix(1,1,values[0]), B, B);
            flip1.insert_block(Matrix(1,1,values[1]), C, C);
            flip1.insert_block(Matrix(1,1,values[2]), B, C);
            flip1.insert_block(Matrix(1,1,values[3]), C, B);

            op_t flip2; flip2.twoS = 2; flip2.twoSaction = -2;
            flip2.insert_block(Matrix(1,1,values[4]), B, B);
            flip2.insert_block(Matrix(1,1,values[5]), C, C);
            flip2.insert_block(Matrix(1,1,values[6]), B, C);
            flip2.insert_block(Matrix(1,1,values[7]), C, B);

            op_t count; count.twoS = 0;
            count.insert_block(Matrix(1,1,2), A, A);
            count.insert_block(Matrix(1,1,1), B, B);
            count.insert_block(Matrix(1,1,1), C, C);

            MPOTensor<Matrix, SymmGroup> op(1,1);
            if (p==i) {
                op = MPOTensor<Matrix, SymmGroup>(1,2);
                op.set(0,0, flip1, -std::sqrt(3.));
                op.set(0,1, count, 0.5);
            }
            else if (p==j) {
                op = MPOTensor<Matrix, SymmGroup>(2,1);
                op.set(0,0, flip2, 1.0);
                op.set(1,0, count, 1.0);
            }
            else if (i < p && p < j)
            {
                op = MPOTensor<Matrix, SymmGroup>(2,2);
                op.set(0,0, fill, 1.0);
                op.set(1,1, ident, 1.0);
            }
            else
                op.set(0,0, ident, 1.0);

            ret[p] = op;
        }
        return ret;
    }

} // namespace SU2

    template<class Matrix>
    MPO<Matrix, TwoU1PG> make_e2d_d2e_abelian(int i, int j, std::vector<int> site_irreps)
    {
        typedef block_matrix<Matrix, TwoU1PG> op_t;
        //boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
        MPO<Matrix, TwoU1PG> ret(site_irreps.size());
        for (int p=0; p<site_irreps.size(); ++p)
        {
            typedef tag_detail::tag_type tag_type;
            typename TwoU1PG::charge A(0), B(0), C(0), D(0);
            B[0] = 1; C[1] = 1;
            D[0] = 1; D[1] = 1;
            B[2] = site_irreps[p];
            C[2] = site_irreps[p];

            block_matrix<Matrix, TwoU1PG> ident;
            ident.insert_block(Matrix(1,1,1), A, A);
            ident.insert_block(Matrix(1,1,1), B, B);
            ident.insert_block(Matrix(1,1,1), C, C);
            ident.insert_block(Matrix(1,1,1), D, D);

            op_t e2d_op;
            e2d_op.insert_block(Matrix(1,1,1), D, A);

            op_t d2e_op;
            d2e_op.insert_block(Matrix(1,1,1), A, D);

            MPOTensor<Matrix, TwoU1PG> op(1,1);
            if (p==i)
                op.set(0,0, e2d_op, 1.0);
            else if (p==j)
                op.set(0,0, d2e_op, 1.0);
            else
                op.set(0,0, ident, 1.0);

            ret[p] = op;
        }
        return ret;
    }

    template<class Matrix>
    MPO<Matrix, TwoU1PG> make_count_abelian(int i, std::vector<int> site_irreps)
    {
        typedef block_matrix<Matrix, TwoU1PG> op_t;
        //boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
        MPO<Matrix, TwoU1PG> ret(site_irreps.size());
        for (int p=0; p<site_irreps.size(); ++p)
        {
            typedef tag_detail::tag_type tag_type;
            typename TwoU1PG::charge A(0), B(0), C(0), D(0);
            B[0] = 1; C[1] = 1;
            D[0] = 1; D[1] = 1;
            B[2] = site_irreps[p];
            C[2] = site_irreps[p];

            block_matrix<Matrix, TwoU1PG> ident;
            ident.insert_block(Matrix(1,1,1), A, A);
            ident.insert_block(Matrix(1,1,1), B, B);
            ident.insert_block(Matrix(1,1,1), C, C);
            ident.insert_block(Matrix(1,1,1), D, D);

            op_t count;
            count.insert_block(Matrix(1,1,1), B, B);
            count.insert_block(Matrix(1,1,1), C, C);
            count.insert_block(Matrix(1,1,2), D, D);

            MPOTensor<Matrix, TwoU1PG> op(1,1);
            if (p==i)
                op.set(0,0, count, 1.0);
            else
                op.set(0,0, ident, 1.0);

            ret[p] = op;
        }
        return ret;
    }

    template<class Matrix>
    MPO<Matrix, TwoU1PG> make_flip_abelian(int i, int j, std::vector<int> site_irreps)
    {
        typedef block_matrix<Matrix, TwoU1PG> op_t;
        //boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
        MPO<Matrix, TwoU1PG> ret(site_irreps.size());
        for (int p=0; p<site_irreps.size(); ++p)
        {
            typedef tag_detail::tag_type tag_type;
            typename TwoU1PG::charge A(0), B(0), C(0), D(0);
            B[0] = 1; C[1] = 1;
            D[0] = 1; D[1] = 1;
            B[2] = site_irreps[p];
            C[2] = site_irreps[p];

            block_matrix<Matrix, TwoU1PG> ident;
            ident.insert_block(Matrix(1,1,1), A, A);
            ident.insert_block(Matrix(1,1,1), B, B);
            ident.insert_block(Matrix(1,1,1), C, C);
            ident.insert_block(Matrix(1,1,1), D, D);

            op_t count_up;
            count_up.insert_block(Matrix(1,1,1), B, B);
            count_up.insert_block(Matrix(1,1,1), D, D);

            op_t count_down;
            count_down.insert_block(Matrix(1,1,1), C, C);
            count_down.insert_block(Matrix(1,1,1), D, D);

            op_t flip_up_down;
            flip_up_down.insert_block(Matrix(1,1,1), B, C);

            op_t flip_down_up;
            flip_down_up.insert_block(Matrix(1,1,1), C, B);

            MPOTensor<Matrix, TwoU1PG> op(1,1);
            if (p==i) {
                op = MPOTensor<Matrix, TwoU1PG>(1,4);
                op.set(0,0, flip_up_down, 1.0);
                op.set(0,1, flip_down_up, 1.0);
                op.set(0,2, count_up, 1.0);
                op.set(0,3, count_down, 1.0);
            }
            else if (p==j) {
                op = MPOTensor<Matrix, TwoU1PG>(4,1);
                op.set(0,0, flip_down_up, 1.0);
                op.set(1,0, flip_up_down, 1.0);
                op.set(2,0, count_up, 1.0);
                op.set(3,0, count_down, 1.0);
            }
            else if (i < p && p < j)
            {
                op = MPOTensor<Matrix, TwoU1PG>(4,4);
                op.set(0,0, ident, 1.0);
                op.set(1,1, ident, 1.0);
                op.set(2,2, ident, 1.0);
                op.set(3,3, ident, 1.0);
            }
            else
                op.set(0,0, ident, 1.0);

            ret[p] = op;
        }
        return ret;
    }

#endif
