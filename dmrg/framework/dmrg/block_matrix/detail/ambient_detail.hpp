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

#define ATOMIC(condition, kernel, ...) assert(condition); ambient::push< ambient::numeric::kernels::kernel ## _atomic<T> >(__VA_ARGS__);

template<class T, class SymmGroup>
class block_matrix;

namespace maquis { namespace dmrg { namespace detail {

    template <typename T>
    void reshape_l2b(ambient::numeric::matrix<T>& out, const ambient::numeric::matrix<T>& in,
                     size_t in_left_offset, size_t in_phys_offset, 
                     size_t out_left_offset, size_t out_right_offset,
                     size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
    {
        ATOMIC(in.atomic() && out.atomic(), reshape_l2b, out, in, in_left_offset, in_phys_offset, 
                                                         out_left_offset, out_right_offset, 
                                                         sdim1, sdim2, ldim, rdim);
    }

    template <typename T>
    void reshape_b2l(ambient::numeric::matrix<T>& out, const ambient::numeric::matrix<T>& in,
                     size_t in_left_offset, size_t in_right_offset, 
                     size_t out_left_offset, size_t out_phys_offset,
                     size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
    {
        ATOMIC(in.atomic() && out.atomic(), reshape_b2l, out, in, in_left_offset, in_right_offset, 
                                                         out_left_offset, out_phys_offset, 
                                                         sdim1, sdim2, ldim, rdim);
    }

    template <typename T>
    inline void reshape_l2r(const ambient::numeric::matrix<T>& left, ambient::numeric::matrix<T>& right,
                            size_t left_offset, size_t right_offset, size_t sdim, size_t ldim, size_t rdim)
    {
        ATOMIC(left.atomic() && right.atomic(), reshape_l2r, left, right, left_offset, right_offset, 
                                                             sdim, ldim, rdim);
    }
    
    template <typename T>
    inline void reshape_r2l(ambient::numeric::matrix<T>& left, const ambient::numeric::matrix<T>& right,
                            size_t left_offset, size_t right_offset, size_t sdim, size_t ldim, size_t rdim)
    {
        ATOMIC(left.atomic() && right.atomic(), reshape_r2l, left, right, left_offset, right_offset, 
                                                             sdim, ldim, rdim);
    }
    
    template <typename T>
    inline void lb_tensor_mpo(ambient::numeric::matrix<T>& out, const ambient::numeric::matrix<T>& in, const ambient::numeric::matrix<T>& alfa,
                              size_t out_offset, size_t in_offset, size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
    {
        ATOMIC(out.atomic() && in.atomic() && alfa.atomic(), lb_tensor_mpo, out, in, alfa, out_offset, in_offset, 
                                                                            sdim1, sdim2, ldim, rdim);
    }
    
    template <typename T>
    inline void rb_tensor_mpo(ambient::numeric::matrix<T>& out, const ambient::numeric::matrix<T>& in, const ambient::numeric::matrix<T>& alfa,
                              size_t out_offset, size_t in_offset, size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
    {
        ATOMIC(out.atomic() && in.atomic() && alfa.atomic(), rb_tensor_mpo, out, in, alfa, out_offset, in_offset, 
                                                                            sdim1, sdim2, ldim, rdim);
    }
   
    template<class T, class SymmGroup>
    std::vector<double> bond_renyi_entropies(const block_matrix<ambient::numeric::diagonal_matrix<T>, SymmGroup>& set){
        std::vector< std::vector<double> > r(set.n_blocks());
        for(std::size_t k = 0; k < set.n_blocks(); ++k){
            std::vector<T>* v_ptr = &r[k];
            ATOMIC(set[k].atomic(), round_square, set[k], v_ptr);
        }
        ambient::playout();

        std::vector<double> ret;
        for(std::size_t k = 0; k < set.n_blocks(); ++k)
            std::copy(r[k].begin(), r[k].end(), std::back_inserter(ret));
        
        return ret;
    }

    template <typename T>
    inline void left_right_boundary_init(ambient::numeric::matrix<T>& a){
        a.fill_value(1.0);
    }
        
} } }

#undef ATOMIC
#endif
