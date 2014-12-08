/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
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
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <sys/stat.h>

#include <boost/lambda/lambda.hpp>

using std::cerr;
using std::cout;
using std::endl;

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>
#include "dmrg/block_matrix/detail/alps.hpp"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/boundary.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/contractions.h"

typedef alps::numeric::matrix<double> matrix;

typedef TwoU1LPG grp;

template<class Matrix>
MPO<Matrix, grp> make_mpo_5(int i, int j, int k, int l, std::vector<int> site_irreps)
{
    typedef block_matrix<Matrix, grp> op_t;
    MPO<Matrix, grp> ret(site_irreps.size());
    for (int p=0; p<site_irreps.size(); ++p)
    {
        typedef tag_detail::tag_type tag_type;
        typename grp::charge A(0), B(0), C(0);
        B[0]=1; C[1]=1;
        B[2] = site_irreps[p];
        C[2] = site_irreps[p];

        block_matrix<Matrix, grp> ident_unbarred;
        ident_unbarred.insert_block(Matrix(1,1,1), A, A);
        ident_unbarred.insert_block(Matrix(1,1,1), B, B);

        block_matrix<Matrix, grp> ident_barred;
        ident_barred.insert_block(Matrix(1,1,1), A, A);
        ident_barred.insert_block(Matrix(1,1,1), C, C);

        block_matrix<Matrix, grp> fill_unbarred;
        fill_unbarred.insert_block(Matrix(1,1,1), A, A);
        fill_unbarred.insert_block(Matrix(1,1,-1), B, B);
 
        block_matrix<Matrix, grp> fill_barred;
        fill_barred.insert_block(Matrix(1,1,1), A, A);
        fill_barred.insert_block(Matrix(1,1,-1), C, C);

        /***********************************************/

        op_t count_unbarred;
        count_unbarred.insert_block(Matrix(1,1,1), B, B);

        op_t count_barred;
        count_barred.insert_block(Matrix(1,1,1), C, C);

        /***********************************************/

        op_t create_unbarred;
        create_unbarred.insert_block(Matrix(1,1,1), A, B);

        op_t create_barred;
        create_barred.insert_block(Matrix(1,1,1), A, C);

        op_t destroy_unbarred;
        destroy_unbarred.insert_block(Matrix(1,1,1), B, A);

        op_t destroy_barred;
        destroy_barred.insert_block(Matrix(1,1,1), C, A);

        /***********************************************/

        op_t destroy_unbarred_fill;
        destroy_unbarred_fill.insert_block(Matrix(1,1,-1), B, A);

        op_t destroy_barred_fill;
        destroy_barred_fill.insert_block(Matrix(1,1,-1), C, A);

        /***********************************************/

        maquis::cout << std::fixed << std::setprecision(10);
        MPOTensor<Matrix, grp> op(1,1);
        if (p==i) {
            op.set(0,0, create_barred, -1.0);
        } else if (p==j) {
            op.set(0,0, destroy_unbarred_fill, 1.0);
        } else if (p==k) {
            op.set(0,0, create_unbarred, 1.0);
        } else if (p==l) {
            op.set(0,0, destroy_barred_fill, 1.0);
        } else if ( (j < p && p < k) || (l < p && p < i) ) {
            if (p < site_irreps.size() / 2)
                op.set(0,0, fill_unbarred, 1.0);
            else
                op.set(0,0, fill_barred, 1.0);
        } else {
            if (p < site_irreps.size() / 2)
                op.set(0,0, ident_unbarred, 1.0);
            else
                op.set(0,0, ident_barred, 1.0);
        }
        ret[p] = op;
    }
    return ret;
}

template<class Matrix>
MPO<Matrix, grp> make_mpo_4(int i, int j, int k, int l, std::vector<int> site_irreps)
{
    typedef block_matrix<Matrix, grp> op_t;
    MPO<Matrix, grp> ret(site_irreps.size());
    for (int p=0; p<site_irreps.size(); ++p)
    {
        typedef tag_detail::tag_type tag_type;
        typename grp::charge A(0), B(0), C(0);
        B[0]=1; C[1]=1;
        B[2] = site_irreps[p];
        C[2] = site_irreps[p];

        block_matrix<Matrix, grp> ident_unbarred;
        ident_unbarred.insert_block(Matrix(1,1,1), A, A);
        ident_unbarred.insert_block(Matrix(1,1,1), B, B);

        block_matrix<Matrix, grp> ident_barred;
        ident_barred.insert_block(Matrix(1,1,1), A, A);
        ident_barred.insert_block(Matrix(1,1,1), C, C);

        block_matrix<Matrix, grp> fill_unbarred;
        fill_unbarred.insert_block(Matrix(1,1,1), A, A);
        fill_unbarred.insert_block(Matrix(1,1,-1), B, B);
 
        block_matrix<Matrix, grp> fill_barred;
        fill_barred.insert_block(Matrix(1,1,1), A, A);
        fill_barred.insert_block(Matrix(1,1,-1), C, C);

        /***********************************************/

        op_t count_unbarred;
        count_unbarred.insert_block(Matrix(1,1,1), B, B);

        op_t count_barred;
        count_barred.insert_block(Matrix(1,1,1), C, C);

        /***********************************************/

        op_t create_unbarred;
        create_unbarred.insert_block(Matrix(1,1,1), A, B);

        op_t create_barred;
        create_barred.insert_block(Matrix(1,1,1), A, C);

        op_t destroy_unbarred;
        destroy_unbarred.insert_block(Matrix(1,1,1), B, A);

        op_t destroy_barred;
        destroy_barred.insert_block(Matrix(1,1,1), C, A);

        /***********************************************/

        op_t destroy_unbarred_fill;
        destroy_unbarred_fill.insert_block(Matrix(1,1,-1), B, A);

        op_t destroy_barred_fill;
        destroy_barred_fill.insert_block(Matrix(1,1,-1), C, A);

        /***********************************************/

        maquis::cout << std::fixed << std::setprecision(10);
        MPOTensor<Matrix, grp> op(1,1);
        if (p==i) {
            op.set(0,0, create_unbarred, 1.0);
        } else if (p==j) {
            op.set(0,0, destroy_barred, 1.0);
        } else if (p==k) {
            op.set(0,0, create_unbarred, 1.0);
        } else if (p==l) {
            op.set(0,0, destroy_barred_fill, -1.0);
        } else if ( (l < p && p < j) || (k < p && p < i) ) {
            if  ( p < site_irreps.size() /2  )
                op.set(0,0, fill_unbarred, 1.0);
            else
                op.set(0,0, fill_barred, 1.0);
        } else {
            if ( p < site_irreps.size() < 2 )
                op.set(0,0, ident_unbarred, 1.0);
            else
                op.set(0,0, ident_barred, 1.0);
        }
        ret[p] = op;
    }
    return ret;
}

template<class Matrix>
MPO<Matrix, grp> make_mpo_3(int i, int j, int k, int l, std::vector<int> site_irreps)
{
    typedef block_matrix<Matrix, grp> op_t;
    MPO<Matrix, grp> ret(site_irreps.size());
    for (int p=0; p<site_irreps.size(); ++p)
    {
        typedef tag_detail::tag_type tag_type;
        typename grp::charge A(0), B(0), C(0);
        B[0]=1; C[1]=1;
        B[2] = site_irreps[p];
        C[2] = site_irreps[p];

        block_matrix<Matrix, grp> ident_unbarred;
        ident_unbarred.insert_block(Matrix(1,1,1), A, A);
        ident_unbarred.insert_block(Matrix(1,1,1), B, B);

        block_matrix<Matrix, grp> ident_barred;
        ident_barred.insert_block(Matrix(1,1,1), A, A);
        ident_barred.insert_block(Matrix(1,1,1), C, C);

        block_matrix<Matrix, grp> fill_unbarred;
        fill_unbarred.insert_block(Matrix(1,1,1), A, A);
        fill_unbarred.insert_block(Matrix(1,1,-1), B, B);
 
        block_matrix<Matrix, grp> fill_barred;
        fill_barred.insert_block(Matrix(1,1,1), A, A);
        fill_barred.insert_block(Matrix(1,1,-1), C, C);

        /***********************************************/

        op_t count_unbarred;
        count_unbarred.insert_block(Matrix(1,1,1), B, B);

        op_t count_barred;
        count_barred.insert_block(Matrix(1,1,1), C, C);

        /***********************************************/

        op_t create_unbarred;
        create_unbarred.insert_block(Matrix(1,1,1), A, B);

        op_t create_barred;
        create_barred.insert_block(Matrix(1,1,1), A, C);

        op_t destroy_unbarred;
        destroy_unbarred.insert_block(Matrix(1,1,1), B, A);

        op_t destroy_barred;
        destroy_barred.insert_block(Matrix(1,1,1), C, A);

        /***********************************************/

        op_t destroy_unbarred_fill;
        destroy_unbarred_fill.insert_block(Matrix(1,1,-1), B, A);

        op_t destroy_barred_fill;
        destroy_barred_fill.insert_block(Matrix(1,1,-1), C, A);

        /***********************************************/

        maquis::cout << std::fixed << std::setprecision(10);
        MPOTensor<Matrix, grp> op(1,1);
        if (p==i) {
            op.set(0,0, create_barred, -1.0);
        } else if (p==j) {
            op.set(0,0, destroy_unbarred, 1.0);
        } else if (p==k) {
            op.set(0,0, create_barred, 1.0);
        } else if (p==l) {
            op.set(0,0, destroy_unbarred_fill, 1.0);
        } else if ( (l < p && p < j) || (k < p && p < i) ) {
            if (p < site_irreps.size() / 2)
                op.set(0,0, fill_unbarred, 1.0);
            else
                op.set(0,0, fill_barred, 1.0);
        } else {
            if (p < site_irreps.size() / 2)
                op.set(0,0, ident_unbarred, 1.0);
            else
                op.set(0,0, ident_barred, 1.0);
        }
        ret[p] = op;
    }
    return ret;
}

template<class Matrix>
MPO<Matrix, grp> make_mpo_2(int i, int j, int k, int l, std::vector<int> site_irreps)
{
    typedef block_matrix<Matrix, grp> op_t;
    MPO<Matrix, grp> ret(site_irreps.size());
    for (int p=0; p<site_irreps.size(); ++p)
    {
        typedef tag_detail::tag_type tag_type;
        typename grp::charge A(0), B(0), C(0);
        B[0]=1; C[1]=1;
        B[2] = site_irreps[p];
        C[2] = site_irreps[p];

        block_matrix<Matrix, grp> ident_unbarred;
        ident_unbarred.insert_block(Matrix(1,1,1), A, A);
        ident_unbarred.insert_block(Matrix(1,1,1), B, B);

        block_matrix<Matrix, grp> ident_barred;
        ident_barred.insert_block(Matrix(1,1,1), A, A);
        ident_barred.insert_block(Matrix(1,1,1), C, C);

        block_matrix<Matrix, grp> fill_unbarred;
        fill_unbarred.insert_block(Matrix(1,1,1), A, A);
        fill_unbarred.insert_block(Matrix(1,1,-1), B, B);
 
        block_matrix<Matrix, grp> fill_barred;
        fill_barred.insert_block(Matrix(1,1,1), A, A);
        fill_barred.insert_block(Matrix(1,1,-1), C, C);

        /***********************************************/

        op_t count_unbarred;
        count_unbarred.insert_block(Matrix(1,1,1), B, B);

        op_t count_barred;
        count_barred.insert_block(Matrix(1,1,1), C, C);

        /***********************************************/

        op_t create_unbarred;
        create_unbarred.insert_block(Matrix(1,1,1), A, B);

        op_t create_barred;
        create_barred.insert_block(Matrix(1,1,1), A, C);

        op_t destroy_unbarred;
        destroy_unbarred.insert_block(Matrix(1,1,1), B, A);

        op_t destroy_barred;
        destroy_barred.insert_block(Matrix(1,1,1), C, A);

        /***********************************************/

        op_t destroy_unbarred_fill;
        destroy_unbarred_fill.insert_block(Matrix(1,1,-1), B, A);

        op_t destroy_barred_fill;
        destroy_barred_fill.insert_block(Matrix(1,1,-1), C, A);

        /***********************************************/

        maquis::cout << std::fixed << std::setprecision(10);
        MPOTensor<Matrix, grp> op(1,1);
        if (p==i) {
            op.set(0,0,count_unbarred, 1.0);
        } else if (p==j) {
            op.set(0,0, destroy_unbarred, 1.0);
        } else if (p==k) {
            op.set(0,0, create_unbarred, -1.0);
        } else if ( (j < p && p < k) || (k < p && p < j) ) {
            if (p < site_irreps.size() / 2)
                op.set(0,0, fill_unbarred, 1.0);
            else
                op.set(0,0, fill_barred, 1.0);
        } else {
            if (p < site_irreps.size() / 2)
                op.set(0,0, ident_unbarred, 1.0);
            else
                op.set(0,0, ident_barred, 1.0);
        }
        ret[p] = op;
    }
    return ret;
}

template<class Matrix>
MPO<Matrix, grp> make_1et(int i, int j, std::vector<int> site_irreps)
{
    typedef block_matrix<Matrix, grp> op_t;
    MPO<Matrix, grp> ret(site_irreps.size());
    for (int p=0; p<site_irreps.size(); ++p)
    {
        typedef tag_detail::tag_type tag_type;
        typename grp::charge A(0), B(0), C(0);
        B[0]=1; C[1]=1;
        B[2] = site_irreps[p];
        C[2] = site_irreps[p];

        block_matrix<Matrix, grp> ident_unbarred;
        ident_unbarred.insert_block(Matrix(1,1,1), A, A);
        ident_unbarred.insert_block(Matrix(1,1,1), B, B);

        block_matrix<Matrix, grp> ident_barred;
        ident_barred.insert_block(Matrix(1,1,1), A, A);
        ident_barred.insert_block(Matrix(1,1,1), C, C);

        block_matrix<Matrix, grp> fill_unbarred;
        fill_unbarred.insert_block(Matrix(1,1,1), A, A);
        fill_unbarred.insert_block(Matrix(1,1,-1), B, B);
 
        block_matrix<Matrix, grp> fill_barred;
        fill_barred.insert_block(Matrix(1,1,1), A, A);
        fill_barred.insert_block(Matrix(1,1,-1), C, C);

        block_matrix<Matrix, grp> transfer;
        transfer.insert_block(Matrix(1,1,1), A, A);
        transfer.insert_block(Matrix(1,1,-1), B, C);

        /***********************************************/

        op_t count_unbarred;
        count_unbarred.insert_block(Matrix(1,1,1), B, B);

        op_t count_barred;
        count_barred.insert_block(Matrix(1,1,1), C, C);

        /***********************************************/

        op_t create_unbarred;
        create_unbarred.insert_block(Matrix(1,1,1), A, B);

        op_t create_barred;
        create_barred.insert_block(Matrix(1,1,1), A, C);

        op_t destroy_unbarred;
        destroy_unbarred.insert_block(Matrix(1,1,1), B, A);

        op_t destroy_barred;
        destroy_barred.insert_block(Matrix(1,1,1), C, A);

        /***********************************************/

        op_t destroy_unbarred_fill;
        destroy_unbarred_fill.insert_block(Matrix(1,1,-1), B, A);

        op_t destroy_barred_fill;
        destroy_barred_fill.insert_block(Matrix(1,1,-1), C, A);

        /***********************************************/

        MPOTensor<Matrix, grp> op(1,1);
        int half_lattice = (site_irreps.size()/2)-1;

        int perm_sign = 1;

        if (p==i) {
            if (i < site_irreps.size() / 2)
                op.set(0,0, create_unbarred, 1.0);
            else
                op.set(0,0, create_barred, 1.0);
        }
        else if (p==j) {
            if (j < i) {
                perm_sign *= -1;
                if (j < site_irreps.size() / 2)
                    op.set(0,0, destroy_unbarred_fill, 1.0*perm_sign);
                else
                    op.set(0,0, destroy_barred_fill, 1.0*perm_sign);
            } else {
                if (j < site_irreps.size() / 2)
                    op.set(0,0, destroy_unbarred, 1.0*perm_sign);
                else
                    op.set(0,0, destroy_barred, 1.0*perm_sign);
            }
        }
        else if ((i < p && p < j) || (j < p && p < i)) {
            if (p < site_irreps.size() / 2)
                op.set(0,0, fill_unbarred, 1.0);
            else
                op.set(0,0, fill_barred, 1.0);
        }     
        else {
            if (p < site_irreps.size() / 2)
                op.set(0,0, ident_unbarred, 1.0);
            else
                op.set(0,0, ident_barred, 1.0);
        }

        ret[p] = op;
    }
    return ret;
}

template<class Matrix>
MPO<Matrix, grp> make_2et(int i, int j, int k, int l,std::vector<int> site_irreps)
{
    typedef block_matrix<Matrix, grp> op_t;
    MPO<Matrix, grp> ret(site_irreps.size());
    for (int p=0; p<site_irreps.size(); ++p)
    {
        typedef tag_detail::tag_type tag_type;
        typename grp::charge A(0), B(0), C(0);
        B[0]=1; C[1]=1;
        B[2] = site_irreps[p];
        C[2] = site_irreps[p];

        block_matrix<Matrix, grp> ident_unbarred;
        ident_unbarred.insert_block(Matrix(1,1,1), A, A);
        ident_unbarred.insert_block(Matrix(1,1,1), B, B);

        block_matrix<Matrix, grp> ident_barred;
        ident_barred.insert_block(Matrix(1,1,1), A, A);
        ident_barred.insert_block(Matrix(1,1,1), C, C);

        block_matrix<Matrix, grp> fill_unbarred;
        fill_unbarred.insert_block(Matrix(1,1,1), A, A);
        fill_unbarred.insert_block(Matrix(1,1,-1), B, B);
 
        block_matrix<Matrix, grp> fill_barred;
        fill_barred.insert_block(Matrix(1,1,1), A, A);
        fill_barred.insert_block(Matrix(1,1,-1), C, C);

        block_matrix<Matrix, grp> transfer;
        transfer.insert_block(Matrix(1,1,1), A, A);
        transfer.insert_block(Matrix(1,1,-1), B, C);

        /***********************************************/

        op_t count_unbarred;
        count_unbarred.insert_block(Matrix(1,1,1), B, B);

        op_t count_barred;
        count_barred.insert_block(Matrix(1,1,1), C, C);

        /***********************************************/

        op_t create_unbarred;
        create_unbarred.insert_block(Matrix(1,1,1), A, B);

        op_t create_barred;
        create_barred.insert_block(Matrix(1,1,1), A, C);

        op_t destroy_unbarred;
        destroy_unbarred.insert_block(Matrix(1,1,1), B, A);

        op_t destroy_barred;
        destroy_barred.insert_block(Matrix(1,1,1), C, A);

        /***********************************************/

        std::vector<int> idx;
        idx.push_back(i);
        idx.push_back(j);
        idx.push_back(k);
        idx.push_back(l);

        int inv_count=0, n=4;
        for(int c1 = 0; c1 < n - 1; c1++)
            for(int c2 = c1+1; c2 < n; c2++)
                if(idx[c1] > idx[c2]) inv_count++;

        std::sort(idx.begin(),idx.end());
        maquis::cout << inv_count << std::endl;
        //if (inv_count % 2)

        MPOTensor<Matrix, grp> op(1,1);
        int half_lattice = (site_irreps.size()/2);
        if (p==idx[0]) {
            if (p < site_irreps.size() / 2)
                op.set(0,0, create_unbarred, 1.0);
            else
                op.set(0,0, create_barred, 1.0);
        }
        else if (p==idx[1]) {
            if (p < site_irreps.size() / 2)
                op.set(0,0, destroy_unbarred, 1.0);
            else
                op.set(0,0, destroy_barred, 1.0);
        }
        else if (p==idx[2]) {
            if (p < site_irreps.size() / 2)
                op.set(0,0, create_unbarred, 1.0);
            else
                op.set(0,0, create_barred, 1.0);
        }
        else if (p==idx[3]) {
            if (p < site_irreps.size() / 2)
                op.set(0,0, destroy_unbarred, 1.0);
            else
                op.set(0,0, destroy_barred, 1.0);
        }
        else if (idx[0] < p && p < idx[1]) {
            if (p < site_irreps.size() / 2)
                op.set(0,0, fill_unbarred, -1.0);
            else
                op.set(0,0, fill_barred, -1.0);
        }     
        else if (idx[2] < p && p < idx[3]) {
            if (p < site_irreps.size() / 2)
                op.set(0,0, fill_unbarred, -1.0);
            else
                op.set(0,0, fill_barred, -1.0);
        }     
        else {
            if (p < site_irreps.size() / 2)
                op.set(0,0, ident_unbarred, 1.0);
            else
                op.set(0,0, ident_barred, 1.0);
        }
        //}
        ret[p] = op;
    }
    return ret;
}
int main(int argc, char ** argv)
{
    try {
        if (argc < 2) {
            std::cout << "Usage: " << argv[0] << " <mps.h5>" << std::endl;
            return 1;
        }
        //cout.precision(5);

        MPS<matrix, grp> mps;
        load(argv[1], mps);
        size_t L = mps.size();

        std::vector<int> site_irreps;
        for (int i=0; i < L; ++i)
            site_irreps.push_back(mps[i].site_dim()[0].first[2]);

        std::cout << "site irreps: ";
        std::copy(site_irreps.begin(), site_irreps.end(), std::ostream_iterator<int>(std::cout, " "));
        std::cout << std::endl;

        //MPO<matrix, grp> mpo = make_1et<matrix>(4,2,site_irreps);
        //double expectation_value = expval(mps, mpo, false);
        //maquis::cout << "5 3: " << expectation_value << std::endl;

        //mpo = make_1et<matrix>(2,4,site_irreps);
        //expectation_value = expval(mps,mpo,false);
        //maquis::cout << "3 5: " << expectation_value << std::endl;
       
        //mpo = make_mpo_2<matrix>(0,2,1,0,site_irreps);
        //expectation_value = expval(mps,mpo,false);
        //maquis::cout << "1 3 2 1: " << expectation_value << std::endl;

        //MPO<matrix, grp> mpo = make_mpo_2<matrix>(0,3,2,0,site_irreps);
        //double expectation_value = expval(mps,mpo,false);
        //maquis::cout << "1 4 3 1: " << expectation_value << std::endl;

        // 4 in 10
        //MPO<matrix,grp> mpo = make_mpo_3<matrix>(7,2,5,0,site_irreps);
        //double expectation_value = expval(mps,mpo,false);
        //maquis::cout << "8 3 6 1: " << expectation_value << std::endl;

        MPO<matrix,grp> mpo = make_mpo_4<matrix>(2,7,0,5,site_irreps);
        double expectation_value = expval(mps,mpo,false);
        maquis::cout << "3 8 1 6: " << expectation_value << std::endl;

        double exp_norm = norm(mps);
        maquis::cout << "norm: " << exp_norm << std::endl;

        //MPO<matrix,grp> mpo = make_mpo_3<matrix>(17,1,15,0,site_irreps);
        //double expectation_value = expval(mps,mpo,false);
        //maquis::cout << "18 2 16 1: " << expectation_value << std::endl;

        //mpo = make_mpo_4<matrix>(6,12,4,9,site_irreps);
        //expectation_value = expval(mps,mpo,false);
        //maquis::cout << "7 13 5 10: " << expectation_value << std::endl;

        //mpo = make_mpo_5<matrix>(13,1,4,9,site_irreps);
        //expectation_value = expval(mps,mpo,false);
        //maquis::cout << "14 2 5 10: " << expectation_value << std::endl;


    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
