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

#ifdef USE_AMBIENT
#include <mpi.h>
#endif
#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#ifdef USE_AMBIENT
    #include "dmrg/block_matrix/detail/ambient.hpp"
    typedef ambient::numeric::tiles<ambient::numeric::matrix<double> > Matrix;
#else
    #include "dmrg/block_matrix/detail/alps.hpp"
    typedef alps::numeric::matrix<double> Matrix;
#endif

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/mp_tensors/twositetensor.h"

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

int main(int argc, char ** argv)
{
    try {
        //if (argc != 2) {
        //    std::cout << "Usage: " << argv[0] << " <mps.h5>" << std::endl;
        //    return 1;
        //}
        //MPS<Matrix, grp> mps;
        //load(argv[1], mps);

        //TwoSiteTensor<Matrix, grp> ts(mps[1], mps[2]);        
        //MPSTensor<Matrix, grp> mpstensor = ts.make_mps();
        //mpstensor.make_left_paired();
        //mpstensor.make_right_paired();
        //ts << mpstensor;

        typename grp::charge A(0), B(0), C(0), D(0);
        A[0] = 2; // 20
        B[0] = 1; B[1] =  1; // 11
        C[0] = 1; C[1] = -1; // 1-1
        // D = 00

        Index<grp> phys;
        phys.insert(std::make_pair(A, 1));
        phys.insert(std::make_pair(B, 1));
        phys.insert(std::make_pair(C, 1));
        phys.insert(std::make_pair(D, 1));

        typedef block_matrix<Matrix, grp> op_t;

        SpinDescriptor<symm_traits::SU2Tag> one_half_up(1,1);
        SpinDescriptor<symm_traits::SU2Tag> one_half_down(1,-1);
        SpinDescriptor<symm_traits::SU2Tag> one_up(2,2);
        SpinDescriptor<symm_traits::SU2Tag> one_flat(2,0);
        SpinDescriptor<symm_traits::SU2Tag> one_down(2,-2);

        // cheaper to use this for spin0 tensors, instead of ident_full
        op_t ident_op;
        ident_op.insert_block(Matrix(1,1,1), A, A);
        ident_op.insert_block(Matrix(1,1,1), B, B);
        ident_op.insert_block(Matrix(1,1,1), C, C);
        ident_op.insert_block(Matrix(1,1,1), D, D);

        // apply if spin > 0
        op_t ident_full_op;
        ident_full_op.insert_block(Matrix(1,1,1), A, A);
        ident_full_op.insert_block(Matrix(1,1,1), D, D);
        ident_full_op.insert_block(Matrix(1,1,1), B, B);
        ident_full_op.insert_block(Matrix(1,1,1), C, C);
        ident_full_op.insert_block(Matrix(1,1,1), B, C);
        ident_full_op.insert_block(Matrix(1,1,1), C, B);

        op_t fill_op;
        fill_op.insert_block(Matrix(1,1,1),  A, A);
        fill_op.insert_block(Matrix(1,1,1),  D, D);
        fill_op.insert_block(Matrix(1,1,-1), B, B);
        fill_op.insert_block(Matrix(1,1,-1), C, C);
        fill_op.insert_block(Matrix(1,1,-1), B, C);
        fill_op.insert_block(Matrix(1,1,-1), C, B);

        /*************************************************************/

        op_t create_fill_op;
        create_fill_op.spin = one_half_up;
        create_fill_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
        create_fill_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);
        create_fill_op.insert_block(Matrix(1,1,1), D, B);
        create_fill_op.insert_block(Matrix(1,1,1), D, C);

        op_t destroy_op;
        destroy_op.spin = one_half_down;
        destroy_op.insert_block(Matrix(1,1,1), A, B);
        destroy_op.insert_block(Matrix(1,1,1), A, C);
        destroy_op.insert_block(Matrix(1,1,sqrt(2.)), B, D);
        destroy_op.insert_block(Matrix(1,1,sqrt(2.)), C, D);

        op_t destroy_fill_op;
        destroy_fill_op.spin = one_half_up;
        destroy_fill_op.insert_block(Matrix(1,1,1), A, B);
        destroy_fill_op.insert_block(Matrix(1,1,1), A, C);
        destroy_fill_op.insert_block(Matrix(1,1,-sqrt(2.)), B, D);
        destroy_fill_op.insert_block(Matrix(1,1,-sqrt(2.)), C, D);

        op_t create_op;
        create_op.spin = one_half_down;
        create_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
        create_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);
        create_op.insert_block(Matrix(1,1,-1), D, B);
        create_op.insert_block(Matrix(1,1,-1), D, C);

        /*************************************************************/

        op_t create_fill_couple_down_op = create_fill_op;
        create_fill_couple_down_op.spin = one_half_down;

        op_t destroy_fill_couple_down_op = destroy_fill_op;
        destroy_fill_couple_down_op.spin = one_half_down;

        op_t create_couple_up_op = create_op;
        create_couple_up_op.spin = one_half_up;

        op_t destroy_couple_up_op = destroy_op;
        destroy_couple_up_op.spin = one_half_up;

        /*************************************************************/

        op_t create_fill_count_op;
        create_fill_count_op.spin = one_half_up;
        create_fill_count_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
        create_fill_count_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);

        op_t destroy_count_op;
        destroy_count_op.spin = one_half_down;
        destroy_count_op.insert_block(Matrix(1,1,1), A, B);
        destroy_count_op.insert_block(Matrix(1,1,1), A, C);

        op_t destroy_fill_count_op;
        destroy_fill_count_op.spin = one_half_up;
        destroy_fill_count_op.insert_block(Matrix(1,1,1), A, B);
        destroy_fill_count_op.insert_block(Matrix(1,1,1), A, C);

        op_t create_count_op;
        create_count_op.spin = one_half_down;
        create_count_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
        create_count_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);

        /*************************************************************/

        op_t count_op;
        count_op.insert_block(Matrix(1,1,2), A, A);
        count_op.insert_block(Matrix(1,1,1), B, B);
        count_op.insert_block(Matrix(1,1,1), C, C);

        op_t docc_op;
        docc_op.insert_block(Matrix(1,1,1), A, A);

        op_t e2d_op;
        e2d_op.insert_block(Matrix(1,1,1), D, A);

        op_t d2e_op;
        d2e_op.insert_block(Matrix(1,1,1), A, D);

        op_t count_fill_op;
        count_fill_op.insert_block(Matrix(1,1,2),  A, A);
        count_fill_op.insert_block(Matrix(1,1,-1), B, B);
        count_fill_op.insert_block(Matrix(1,1,-1), C, C);
        count_fill_op.insert_block(Matrix(1,1,-1), B, C);
        count_fill_op.insert_block(Matrix(1,1,-1), C, B);

        op_t flip_to_S2_op;
        flip_to_S2_op.spin = one_up;
        flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2)), B, B);
        flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2.)), C, C);
        flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2.)),  B, C);
        flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2.)),  C, B);

        op_t flip_to_S0_op = flip_to_S2_op;
        flip_to_S0_op.spin = one_down;

        op_t flip_S0_op = flip_to_S2_op;
        flip_S0_op.spin = one_flat;

        /*************************************************************/

        op_t kronop;

        op_kron(phys, phys, count_op, count_op, kronop) ;
        maquis::cout << "Kronecker product count * count\n";
        maquis::cout << kronop << std::endl;

        op_kron(phys, phys, ident_op, create_fill_op, kronop) ;
        maquis::cout << "Kronecker product I * create_fill\n";
        maquis::cout << kronop << std::endl;

        op_kron(phys, phys, ident_op, destroy_fill_op, kronop) ;
        maquis::cout << "Kronecker product I * destroy_fill\n";
        maquis::cout << kronop << std::endl;

        op_kron(phys, phys, create_op, ident_op, kronop) ;
        maquis::cout << "Kronecker product create * I\n";
        maquis::cout << kronop << std::endl;

        op_kron(phys, phys, destroy_op, ident_op, kronop) ;
        maquis::cout << "Kronecker product destroy * I\n";
        maquis::cout << kronop << std::endl;

        op_kron(phys, phys, create_fill_op, destroy_op, kronop) ;
        maquis::cout << "Kronecker product create_fill * destroy\n";
        maquis::cout << kronop << std::endl;

        op_kron(phys, phys, destroy_fill_op, create_op, kronop) ;
        maquis::cout << "Kronecker product destroy_fill * create\n";
        maquis::cout << kronop << std::endl;
            
    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
