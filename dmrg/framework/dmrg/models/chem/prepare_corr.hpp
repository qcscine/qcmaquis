#ifndef PREPARE_OPS_HPP
#define PREPARE_OPS_HPP

#include <vector>

#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/block_matrix/block_matrix.h"

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
            fill.insert_block(Matrix(1,1,1), B, C);  // -1
            fill.insert_block(Matrix(1,1,1), C, B);  // -1

            block_matrix<Matrix, SymmGroup> create;
            create.insert_block(Matrix(1,1,-2.), B, A);       // 1
            create.insert_block(Matrix(1,1,2.), C, A);       // 1
            create.insert_block(Matrix(1,1,-sqrt(2.)), D, B); // 1
            create.insert_block(Matrix(1,1,sqrt(2.)), D, C); // 1

            block_matrix<Matrix, SymmGroup> destroy;
            destroy.insert_block(Matrix(1,1,1), A, B);        // 1 
            destroy.insert_block(Matrix(1,1,-1), A, C);        // 1
            destroy.insert_block(Matrix(1,1,sqrt(2.)), B, D); // 1
            destroy.insert_block(Matrix(1,1,-sqrt(2.)), C, D); // 1

            //tag_type ident = tag_handler->register_op(identity, tag_detail::bosonic);
            MPOTensor<Matrix, SymmGroup> op(1,1);

            if (p == i) {
                //block_matrix<Matrix, SymmGroup> tmp;
                //SU2::gemm(fill, create, tmp);
                //maquis::cout << tmp << std::endl;
                op.set(0,0,destroy,1.0);
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

            block_matrix<Matrix, SymmGroup> create1;
            create1.insert_block(Matrix(1,1,sqrt(2.)), B, A);      
            create1.insert_block(Matrix(1,1,sqrt(2.)), C, A);      
            create1.insert_block(Matrix(1,1,1), D, B);
            create1.insert_block(Matrix(1,1,1), D, C);

            block_matrix<Matrix, SymmGroup> create2;
            create2.insert_block(Matrix(1,1,sqrt(2.)), B, A);      
            create2.insert_block(Matrix(1,1,sqrt(2.)), C, A);      
            create2.insert_block(Matrix(1,1,1), D, B);
            create2.insert_block(Matrix(1,1,1), D, C);

            block_matrix<Matrix, SymmGroup> destroy1;
            destroy1.insert_block(Matrix(1,1,1), A, B);      
            destroy1.insert_block(Matrix(1,1,-1), A, C);     
            destroy1.insert_block(Matrix(1,1,sqrt(2.)), B, D); 
            destroy1.insert_block(Matrix(1,1,-sqrt(2.)), C, D);


            block_matrix<Matrix, SymmGroup> destroy2;
            destroy2.insert_block(Matrix(1,1,1), A, B);        
            destroy2.insert_block(Matrix(1,1,1), A, C);       
            destroy2.insert_block(Matrix(1,1,sqrt(2.)), B, D); 
            destroy2.insert_block(Matrix(1,1,sqrt(2.)), C, D);

            //tag_type ident = tag_handler->register_op(identity, tag_detail::bosonic);
            MPOTensor<Matrix, SymmGroup> op(1,1);

            if (p == i)
                op.set(0,0,create1, 1.0);
            else if (p == j)
                op.set(0,0,create2, 1.0);
            else if (p == k)
                op.set(0,0,destroy1, 1.0);
            else if (p == l)
                op.set(0,0,destroy2, 1.0);
            //else if ( i < p && p < j)
            //    op.set(0,0,fill, 1.0);
            else 
                op.set(0,0,identity, 1.0);

            ret[p] = op;
        }
        return ret;
    }

    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup> make_2rdm_term_custom(int i, int j, int k, int l, int m1, int m2, int m3, int m4, std::vector<int> site_irreps)
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

            block_matrix<Matrix, SymmGroup> create1;
            create1.insert_block(Matrix(1,1,sqrt(2.)), B, A);      
            create1.insert_block(Matrix(1,1,sqrt(2.)), C, A);      
            create1.insert_block(Matrix(1,1,1), D, B);
            create1.insert_block(Matrix(1,1,1), D, C);

            block_matrix<Matrix, SymmGroup> create2;
            create2.insert_block(Matrix(1,1,sqrt(2.)), B, A);      
            create2.insert_block(Matrix(1,1,sqrt(2.)), C, A);      
            create2.insert_block(Matrix(1,1,1), D, B);
            create2.insert_block(Matrix(1,1,1), D, C);

            block_matrix<Matrix, SymmGroup> destroy1;
            destroy1.insert_block(Matrix(1,1,1), A, B);      
            destroy1.insert_block(Matrix(1,1,1), A, C);     
            destroy1.insert_block(Matrix(1,1,sqrt(2.)), B, D); 
            destroy1.insert_block(Matrix(1,1,sqrt(2.)), C, D);


            block_matrix<Matrix, SymmGroup> destroy2;
            destroy2.insert_block(Matrix(1,1,1), A, B);        
            destroy2.insert_block(Matrix(1,1,1), A, C);       
            destroy2.insert_block(Matrix(1,1,sqrt(2.)), B, D); 
            destroy2.insert_block(Matrix(1,1,sqrt(2.)), C, D);


            block_matrix<Matrix, SymmGroup> create1m;
            create1m.insert_block(Matrix(1,1,sqrt(2.)), B, A);      
            create1m.insert_block(Matrix(1,1,-sqrt(2.)), C, A);      
            create1m.insert_block(Matrix(1,1,1), D, B);
            create1m.insert_block(Matrix(1,1,-1), D, C);

            block_matrix<Matrix, SymmGroup> create2m;
            create2m.insert_block(Matrix(1,1,sqrt(2.)), B, A);      
            create2m.insert_block(Matrix(1,1,-sqrt(2.)), C, A);      
            create2m.insert_block(Matrix(1,1,1), D, B);
            create2m.insert_block(Matrix(1,1,-1), D, C);

            block_matrix<Matrix, SymmGroup> destroy1m;
            destroy1m.insert_block(Matrix(1,1,1), A, B);      
            destroy1m.insert_block(Matrix(1,1,-1), A, C);     
            destroy1m.insert_block(Matrix(1,1,sqrt(2.)), B, D); 
            destroy1m.insert_block(Matrix(1,1,-sqrt(2.)), C, D);


            block_matrix<Matrix, SymmGroup> destroy2m;
            destroy2m.insert_block(Matrix(1,1,1), A, B);        
            destroy2m.insert_block(Matrix(1,1,-1), A, C);       
            destroy2m.insert_block(Matrix(1,1,sqrt(2.)), B, D); 
            destroy2m.insert_block(Matrix(1,1,-sqrt(2.)), C, D);

            //tag_type ident = tag_handler->register_op(identity, tag_detail::bosonic);
            MPOTensor<Matrix, SymmGroup> op(1,1);

            if (p == i)
                op.set(0,0,(m1) ? create1m : create1, 1.0);
            else if (p == j)
                op.set(0,0,(m2) ? create2m : create2, 1.0);
            else if (p == k)
                op.set(0,0,(m3) ? destroy1m : destroy1, 1.0);
            else if (p == l)
                op.set(0,0,(m4) ? destroy2m : destroy2, 1.0);
            //else if ( i < p && p < j)
            //    op.set(0,0,fill, 1.0);
            else 
                op.set(0,0,identity, 1.0);

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

} // namespace SU2

#endif
