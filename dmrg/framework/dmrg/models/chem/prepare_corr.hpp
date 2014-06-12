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
    MPO<Matrix, SymmGroup> make_op(int i, int j, std::vector<int> site_irreps)
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
            fill.insert_block(Matrix(1,1,-1), B, B);
            fill.insert_block(Matrix(1,1,-1), C, C);
            fill.insert_block(Matrix(1,1,1), D, D);

            block_matrix<Matrix, SymmGroup> create;
            create.insert_block(Matrix(1,1,1), B, A);
            create.insert_block(Matrix(1,1,1), C, A);
            create.insert_block(Matrix(1,1,1), D, B);
            create.insert_block(Matrix(1,1,1), D, C);

            block_matrix<Matrix, SymmGroup> destroy;
            destroy.insert_block(Matrix(1,1,1), A, B);
            destroy.insert_block(Matrix(1,1,1), A, C);
            destroy.insert_block(Matrix(1,1,1), B, D);
            destroy.insert_block(Matrix(1,1,1), C, D);

            //tag_type ident = tag_handler->register_op(identity, tag_detail::bosonic);
            MPOTensor<Matrix, SymmGroup> op(1,1);

            if (p == i)
                op.set(0,0,destroy);
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
    MPO<Matrix, SymmGroup> make_op_fill(std::vector<int> site_irreps)
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

            block_matrix<Matrix, SymmGroup> fill;
            fill.insert_block(Matrix(1,1,1), A, A);
            fill.insert_block(Matrix(1,1,-1), B, B);
            fill.insert_block(Matrix(1,1,-1), C, C);
            fill.insert_block(Matrix(1,1,1), D, D);
            MPOTensor<Matrix, SymmGroup> op(1,1);

            op.set(0,0,fill, 1.0);

            ret[p] = op;
        }
        return ret;
    }

} // namespace SU2

#endif