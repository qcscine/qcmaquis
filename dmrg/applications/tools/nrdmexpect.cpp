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

typedef TwoU1PG grp;

    template<class Matrix>
    MPO<Matrix, TwoU1PG> make_2rdm_term(int i, int j, int k, int l, std::vector<int> site_irreps)
    {
        typedef block_matrix<Matrix, TwoU1PG> op_t;
        MPO<Matrix, TwoU1PG> ret(site_irreps.size());
        for (int p=0; p<site_irreps.size(); ++p)
        {
            typedef tag_detail::tag_type tag_type;
            typename TwoU1PG::charge A(0), B(0), C(0), D(0);
            A[0] = 1; A[1] = 1;
            B[0] = 1; C[1] = 1;
            B[2] = site_irreps[p];
            C[2] = site_irreps[p];

            block_matrix<Matrix, TwoU1PG> ident;
            ident.insert_block(Matrix(1,1,1), A, A);
            ident.insert_block(Matrix(1,1,1), B, B);
            ident.insert_block(Matrix(1,1,1), C, C);
            ident.insert_block(Matrix(1,1,1), D, D);

            block_matrix<Matrix, TwoU1PG> fill;
            fill.insert_block(Matrix(1,1,1), A, A);
            fill.insert_block(Matrix(1,1,-1), B, B);
            fill.insert_block(Matrix(1,1,-1), C, C);
            fill.insert_block(Matrix(1,1,1), D, D);

            op_t create_up;
            create_up.insert_block(Matrix(1,1,1), C, A);
            create_up.insert_block(Matrix(1,1,1), D, B);
            op_t create_up_fill;
            create_up_fill.insert_block(Matrix(1,1,-1), C, A);
            create_up_fill.insert_block(Matrix(1,1,1), D, B);

            op_t create_down;
            create_down.insert_block(Matrix(1,1,-1), B, A);
            create_down.insert_block(Matrix(1,1,1), D, C);
            op_t create_down_fill;
            create_down_fill.insert_block(Matrix(1,1,1), B, A);
            create_down_fill.insert_block(Matrix(1,1,1), D, C);

            op_t destroy_up;
            destroy_up.insert_block(Matrix(1,1,1), A, C);
            destroy_up.insert_block(Matrix(1,1,1), B, D);
            op_t destroy_up_fill;
            destroy_up_fill.insert_block(Matrix(1,1,1), A, C);
            destroy_up_fill.insert_block(Matrix(1,1,-1), B, D);

            op_t destroy_down;
            destroy_down.insert_block(Matrix(1,1,-1), A, B);
            destroy_down.insert_block(Matrix(1,1,1), C, D);
            op_t destroy_down_fill;
            destroy_down_fill.insert_block(Matrix(1,1,-1), A, B);
            destroy_down_fill.insert_block(Matrix(1,1,-1), C, D);

            MPOTensor<Matrix, TwoU1PG> op(1,1);
            if (p==i) {
                op = MPOTensor<Matrix, TwoU1PG>(1,4);
                op.set(0,0, create_up_fill, 1.0);
                op.set(0,1, create_up_fill, 1.0);
                op.set(0,2, create_down_fill, 1.0);
                op.set(0,3, create_down_fill, 1.0);
            }
            else if (p==j) {
                op = MPOTensor<Matrix, TwoU1PG>(4,4);
                op.set(0,0, create_up, 1.0);
                op.set(1,1, create_down, 1.0);
                op.set(2,2, create_up, 1.0);
                op.set(3,3, create_down, 1.0);
            }
            else if (p==k) {
                op = MPOTensor<Matrix, TwoU1PG>(4,4);
                op.set(0,0, destroy_up_fill, 1.0);
                op.set(1,1, destroy_down_fill, 1.0);
                op.set(2,2, destroy_up_fill, 1.0);
                op.set(3,3, destroy_down_fill, 1.0);
            }
            else if (p==l) {
                op = MPOTensor<Matrix, TwoU1PG>(4,1);
                op.set(0,0, destroy_up, 1.0);
                op.set(1,0, destroy_up, 1.0);
                op.set(2,0, destroy_down, 1.0);
                op.set(3,0, destroy_down, 1.0);
            }
            else if ((i < p && p < j) || (k < p && p < l)) {
                // position is in between first or second pair of operators -> push fill
                op.set(0,0, fill, 1.0);
            }
            else
                op.set(0,0, ident, 1.0);

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
            site_irreps.push_back(mps[i].site_dim()[1].first[2]);

        std::cout << "site irreps: ";
        std::copy(site_irreps.begin(), site_irreps.end(), std::ostream_iterator<int>(std::cout, " "));
        std::cout << std::endl;

        MPO<matrix, grp> mpo = make_2rdm_term<matrix>(0,0,0,0, site_irreps);
        double expectation_value = expval(mps, mpo);
        maquis::cout << "expectation value: " << expectation_value << std::endl;

    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
