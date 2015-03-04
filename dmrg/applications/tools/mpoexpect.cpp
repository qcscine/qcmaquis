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

#include <iostream>
#include <vector>

#include <alps/numeric/matrix.hpp>
#include "dmrg/block_matrix/detail/alps.hpp"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

typedef alps::numeric::matrix<double> matrix;
typedef U1DG grp;

template<class Matrix>
MPO<Matrix, grp> make_1et(int i, int j, std::vector<int> site_irreps)
{
    typedef block_matrix<Matrix, grp> op_t;
    MPO<Matrix, grp> ret(site_irreps.size());
	std::vector<int> op_string;
	op_string.push_back(i); op_string.push_back(j);
	std::sort(op_string.begin(), op_string.end());
	for (int p=0; p<site_irreps.size(); ++p)
    {
        typedef tag_detail::tag_type tag_type;
        typename grp::charge A(0), B(0);
        B[0]=1;
        B[1] = site_irreps[p];

        op_t ident;
        ident.insert_block(Matrix(1,1,1), A, A);
        ident.insert_block(Matrix(1,1,1), B, B);

        op_t fill;
        fill.insert_block(Matrix(1,1,1), A, A);
        fill.insert_block(Matrix(1,1,-1), B, B);

        /***********************************************/

        op_t create;
        create.insert_block(Matrix(1,1,1), A, B);

        op_t destroy;
        destroy.insert_block(Matrix(1,1,1), B, A);

		op_t count;
		count.insert_block(Matrix(1,1,1), B, B);

        /***********************************************/

        op_t destroy_fill;
        destroy_fill.insert_block(Matrix(1,1,-1), B, A);

        /***********************************************/

        maquis::cout << std::fixed << std::setprecision(10);
        MPOTensor<Matrix, grp> op(1,1);
		if ( (i == j) && (p == i) )
			op.set(0,0, count, 1.0);
		else if ( (i != j) && (p == i) )
            op.set(0,0, create, 1.0);
        else if ( (i != j) && (p == j) )
            op.set(0,0, destroy, 1.0);
        else if ( (op_string[0] < p && p < op_string[1]) )
            op.set(0,0, fill, 1.0);
        else
            op.set(0,0, ident, 1.0);
        
        ret[p] = op;
    }
    return ret;
}

template<class Matrix>
MPO<Matrix, grp> make_2et(int i, int j, int k, int l, std::vector<int> site_irreps)
{
    typedef block_matrix<Matrix, grp> op_t;
    MPO<Matrix, grp> ret(site_irreps.size());
	std::vector<int> op_string;
	op_string.push_back(i); op_string.push_back(j); op_string.push_back(k); op_string.push_back(l);
	std::sort(op_string.begin(), op_string.end());
	for (int p=0; p<site_irreps.size(); ++p)
    {
        typedef tag_detail::tag_type tag_type;
        typename grp::charge A(0), B(0);
        B[0]=1;
        B[1] = site_irreps[p];

        op_t ident;
        ident.insert_block(Matrix(1,1,1), A, A);
        ident.insert_block(Matrix(1,1,1), B, B);

        op_t fill;
        fill.insert_block(Matrix(1,1,1), A, A);
        fill.insert_block(Matrix(1,1,-1), B, B);

        /***********************************************/

        op_t create;
        create.insert_block(Matrix(1,1,1), A, B);

        op_t destroy;
        destroy.insert_block(Matrix(1,1,1), B, A);

		op_t count;
		count.insert_block(Matrix(1,1,1), B, B);

        /***********************************************/

        op_t destroy_fill;
        destroy_fill.insert_block(Matrix(1,1,-1), B, A);

        /***********************************************/

        maquis::cout << std::fixed << std::setprecision(10);
        MPOTensor<Matrix, grp> op(1,1);
        if (p==i) {
            op.set(0,0, create, 1.0);
        } else if (p==j) {
            op.set(0,0, destroy, 1.0);
        } else if (p==k) {
            op.set(0,0, create, 1.0);
        } else if (p==l) {
            op.set(0,0, destroy, 1.0);
        } else if ( (j < p && p < l) || (k < p && p < i) ) {
            op.set(0,0, fill, 1.0);
        } else {
            op.set(0,0, ident, 1.0);
        }
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
            site_irreps.push_back(mps[i].site_dim()[0].first[1]);

		// 2 RDM
		//int i,j,k,l;
		//maquis::cout << "Enter 2-RDM element: ";
		//std::cin >> i >> j >> k >> l;
        //MPO<matrix, grp> mpo = make_2et<matrix>(i-1,j-1,k-1,l-1,site_irreps);
        //double expectation_value = expval(mps, mpo, false);
        //maquis::cout << "expectation value: " << expectation_value << std::endl;

		// 1 RDM
		int i,j;
		maquis::cout << "Enter 1-RDM element: ";
		std::cin >> i >> j;
        MPO<matrix, grp> mpo = make_1et<matrix>(i-1,j-1,site_irreps);
        double expectation_value = expval(mps, mpo, false);
        maquis::cout << "expectation value: " << expectation_value << std::endl;

        double exp_norm = norm(mps);
        maquis::cout << "norm: " << exp_norm << std::endl;

    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
