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
MPO<Matrix, grp> make_mpo_2(std::vector<int> site_irreps)
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

        //MPOTensor<Matrix, grp> op(1,1);
        if (p==0) {
            MPOTensor<Matrix, grp> op(1,2);
            op.set(0,0,create_unbarred, 1.0);
            op.set(0,1,destroy_unbarred, 1.0);
        ret[p] = op;
        } else if (p==1) {
            MPOTensor<Matrix, grp> op(2,1);
            op.set(0,0, destroy_unbarred, 1.0);
            op.set(1,0, create_unbarred, 1.0);
        ret[p] = op;
        } else {
            MPOTensor<Matrix, grp> op(1,1);
            op.set(0,0, ident_barred, 1.0);
        ret[p] = op;
        }

    }
    return ret;
}

template<class Matrix>
MPO<Matrix, grp> make_mpo(int i, int j, std::vector<int> site_irreps)
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

        op_t create_unbarred_fill;
        create_unbarred_fill.insert_block(Matrix(1,1,1), A, B);

        op_t create_barred_fill;
        create_barred_fill.insert_block(Matrix(1,1,1), A, C);

        op_t destroy_unbarred_fill;
        destroy_unbarred_fill.insert_block(Matrix(1,1,1), B, A);

        op_t destroy_barred_fill;
        destroy_barred_fill.insert_block(Matrix(1,1,1), C, A);

        /***********************************************/

        MPOTensor<Matrix, grp> op(1,1);
        int half_lattice = (site_irreps.size()/2)-1;
        //if (p==half_lattice) {
        //    if (p==i)
        //        op.set(0,0, gemm(create_unbarred,transfer), 1.0);
        //}
        //else {
        if (p==i) {
            if (i < site_irreps.size() / 2)
                op.set(0,0, create_unbarred, 1.0);
            else
                op.set(0,0, create_barred, 1.0);
        }
        else if (p==j) {
            if (j < site_irreps.size() / 2)
                op.set(0,0, destroy_unbarred, 1.0);
            else
                op.set(0,0, destroy_barred, 1.0);
        }
        else if (i < p && p < j) {
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

        //int pos1=0;
        //int pos2=1;
        //MPO<matrix, grp> mpo = make_mpo<matrix>(pos1,pos2, site_irreps);
        //double expectation_value = expval(mps, mpo, true);
        //maquis::cout << "expectation value for (" << pos1 << "," << pos2 << "): " << expectation_value << std::endl;

        /*
        MPO<matrix, grp> mpo = make_mpo_2<matrix>(site_irreps);
        double expectation_value = expval(mps, mpo, true);
        maquis::cout << "expectation value: " << expectation_value << std::endl;
        */
        MPO<matrix, grp> mpo = make_mpo<matrix>(1,2,site_irreps);
        double expectation_value = expval(mps, mpo, true);
        maquis::cout << "expectation value (2,3): " << expectation_value << std::endl;
        
        mpo = make_mpo<matrix>(0,2,site_irreps);
        expectation_value = expval(mps, mpo, true);
        maquis::cout << "expectation value (1,3): " << expectation_value << std::endl;
        /*
        for (int i=0; i < L; ++i) {
            for (int j=i+1; j < L; ++j) {
                MPO<matrix, grp> mpo = make_mpo<matrix>(i,j, site_irreps);
                double expectation_value = expval(mps, mpo, true);
                maquis::cout << "expectation value for (" << i+1 << "," << j+1 << "): " << expectation_value << std::endl;
            }
        }
        */

    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
