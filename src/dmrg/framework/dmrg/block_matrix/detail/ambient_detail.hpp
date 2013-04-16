/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2012 by Alexandr Kosenkov <alex.kosenkov@gmail.com>
 *                            Timothee Ewart <timothee.ewart@gmail.com>
 *
 *****************************************************************************/

#ifndef MAQUIS_BLOCK_MATRIX_DEATAIL_AMBIENT_MATRIX_DETAIL_HPP
#define MAQUIS_BLOCK_MATRIX_DEATAIL_AMBIENT_MATRIX_DETAIL_HPP

#define value_type typename M::value_type

template<class T, class SymmGroup>
class block_matrix;

namespace maquis { namespace dmrg { namespace detail {

    using ambient::numeric::tiles;
    using ambient::complexity;

    template <typename M>
    inline void op_kron(ambient::numeric::tiles<M>& out, const ambient::numeric::tiles<M>& in, const ambient::numeric::tiles<M>& alfa,
                        size_t out_y_offset, size_t out_x_offset, 
                        size_t ldim1, size_t ldim2, 
                        size_t rdim1, size_t rdim2)
    {
        for(size_t l1 = 0; l1 < ldim1; ++l1)
        for(size_t r1 = 0; r1 < rdim1; ++r1)
        copy_block_s(in, 0, 0, out, out_y_offset + l1*ldim2, out_x_offset + r1*rdim2, 
                     alfa, l1, r1, ldim2, rdim2);
    }

    template <class M>
    inline void reshape_l2b(ambient::numeric::tiles<M>& out, const ambient::numeric::tiles<M>& in,
                     size_t in_left_offset, size_t in_phys_offset, 
                     size_t out_left_offset, size_t out_right_offset,
                     size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
    {
        size_t in_y_offset  = in_left_offset + ldim*in_phys_offset;
        size_t out_y_offset = out_left_offset;

        for(size_t ss1 = 0; ss1 < sdim1; ++ss1){
            for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
                copy_block(in,  in_y_offset, 0, 
                           out, out_y_offset, out_right_offset + rdim*ss2,
                           ldim, rdim);
                in_y_offset += ldim;
            }
            out_y_offset += ldim;
        }
    }

    template <class M>
    inline void reshape_b2l(ambient::numeric::tiles<M>& out, const ambient::numeric::tiles<M>& in,
                     size_t in_left_offset, size_t in_right_offset, 
                     size_t out_left_offset, size_t out_phys_offset,
                     size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
    {
        size_t in_y_offset  = in_left_offset;
        size_t out_y_offset = out_left_offset + out_phys_offset*ldim;

        for(size_t ss1 = 0; ss1 < sdim1; ++ss1){
            for(size_t ss2 = 0; ss2 < sdim2; ++ss2)
            {
                copy_block(in, in_y_offset, in_right_offset + rdim*ss2,
                           out, out_y_offset, 0, 
                           ldim, rdim);
                out_y_offset += ldim;
            }
            in_y_offset += ldim;
        }
    }

    template <class M>
    inline void reshape_l2r(const ambient::numeric::tiles<M>& left, ambient::numeric::tiles<M>& right,
                            size_t left_offset, size_t right_offset, size_t sdim, size_t ldim, size_t rdim)
    {
        for(size_t ss = 0; ss < sdim; ++ss){
        copy_block(left, ss*ldim + left_offset, 0, 
                   right, 0, ss*rdim + right_offset, 
                   ldim, rdim);
        }
    }
    
    template <class M>
    inline void reshape_r2l(ambient::numeric::tiles<M>& left, const ambient::numeric::tiles<M>& right,
                            size_t left_offset, size_t right_offset, size_t sdim, size_t ldim, size_t rdim)
    {
        for(size_t ss = 0; ss < sdim; ++ss)
        copy_block(right, 0, ss*rdim + right_offset, 
                   left, ss*ldim + left_offset, 0, 
                   ldim, rdim);
    }
    
    template <class M1, class M2, class M3>
    inline void lb_tensor_mpo(ambient::numeric::tiles<M1>& out, const ambient::numeric::tiles<M2>& in, const ambient::numeric::tiles<M3>& alfa,
                              size_t out_offset, size_t in_offset, size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
    {
        for(size_t ss2 = 0; ss2 < sdim2; ++ss2)
        for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
        copy_block_sa(in, in_offset + ss1*ldim, 0,
                      out, out_offset + ss2*ldim, 0,
                      alfa, ss1, ss2, ldim, rdim);
    }
    
    template <class M1, class M2, class M3>
    inline void rb_tensor_mpo(ambient::numeric::tiles<M1>& out, const ambient::numeric::tiles<M2>& in, const ambient::numeric::tiles<M3>& alfa,
                              size_t out_offset, size_t in_offset, size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
    {
        for(size_t ss2 = 0; ss2 < sdim2; ++ss2)
        for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
        copy_block_sa(in, 0, in_offset + ss1*rdim,
                      out, 0, out_offset + ss2*rdim,
                      alfa, ss1, ss2, ldim, rdim);
    }
   
    template<class M, class SymmGroup>
    std::vector<double> bond_renyi_entropies(const block_matrix<ambient::numeric::tiles<M>, SymmGroup>& set){
        size_t nblocks = 0;
        size_t k_max = set.n_blocks();
        for(size_t k = 0; k < k_max; ++k)
            nblocks += set[k].data.size();

        size_t vi = 0;
        std::vector< std::vector<double> > r(nblocks);
        for(size_t k = 0; k < k_max; ++k){
            for(size_t kk = 0; kk < set[k].data.size(); kk++){
                std::vector<value_type>* v_ptr = &r[vi++];
                ambient::numeric::kernels::round_square<value_type>::spawn<complexity::N>(set[k][kk], v_ptr);
            }
        }
        ambient::sync();

        std::vector<double> ret;
        for(size_t k = 0; k < nblocks; ++k)
            std::copy(r[k].begin(), r[k].end(), std::back_inserter(ret));
        
        return ret;
    }

    template <class M>
    inline void left_right_boundary_init(ambient::numeric::tiles<M>& a){
        size_t k_max = a.data.size();
        for(size_t k = 0; k < k_max; ++k)
            fill_value(a[k], value_type(1.0));
    }
        
} } }

#undef value_type
#endif
