/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2012 by Alexandr Kosenkov <alex.kosenkov@gmail.com>
 *                            Timothee Ewart <timothee.ewart@gmail.com>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_DETAIL_AMBIENT_HPP
#define MAQUIS_DMRG_DETAIL_AMBIENT_HPP

#include <ambient/numeric/matrix.hpp>
#define ATOMIC(condition, kernel, ...) assert(condition); ambient::push< ambient::numeric::kernels::kernel ## _atomic<T> >(__VA_ARGS__);

namespace maquis { namespace dmrg { namespace detail {
        
        template <typename T>
        inline void reshape_l2r(const ambient::numeric::matrix<T>& left, ambient::numeric::matrix<T>& right,
                                size_t left_offset, size_t right_offset, 
                                size_t sdim, size_t ldim, size_t rdim)
        {
            ATOMIC(left.atomic() && right.atomic(), reshape_l2r, left, right, left_offset, right_offset, sdim, ldim, rdim);
        }
        
        template <typename T>
        inline void reshape_r2l(ambient::numeric::matrix<T>& left, const ambient::numeric::matrix<T>& right,
                                size_t left_offset, size_t right_offset, 
                                size_t sdim, size_t ldim, size_t rdim)
        {
            ATOMIC(left.atomic() && right.atomic(), reshape_r2l, left, right, left_offset, right_offset, sdim, ldim, rdim);
        }
        
        template <typename T>
        inline void lb_tensor_mpo(ambient::numeric::matrix<T>& out, const ambient::numeric::matrix<T>& in, const ambient::numeric::matrix<T>& alfa,
                                  size_t out_offset, size_t in_offset, 
                                  size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
        {
            ATOMIC(out.atomic() && in.atomic() && alfa.atomic(), lb_tensor_mpo, out, in, alfa, out_offset, in_offset, sdim1, sdim2, ldim, rdim);
        }
        
        template <typename T>
        inline void rb_tensor_mpo(ambient::numeric::matrix<T>& out, const ambient::numeric::matrix<T>& in, const ambient::numeric::matrix<T>& alfa,
                                  size_t out_offset, size_t in_offset, 
                                  size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
        {
            ATOMIC(out.atomic() && in.atomic() && alfa.atomic(), rb_tensor_mpo, 
                   out, in, alfa, out_offset, in_offset, sdim1, sdim2, ldim, rdim);
        }
        
        template <typename T>
        inline void bond_renyi_entropies(const ambient::numeric::diagonal_matrix<T>& m, typename alps::numeric::associated_real_vector<ambient::numeric::matrix<T> >::type& sv){
            std::vector<T>* sc_ptr = &sv;
            ambient::push< ambient::numeric::kernels::push_back_sqr_gt<T> >(m, sc_ptr);
            ambient::playout();
        }
        
        template <typename T>
        inline void left_right_boundary_init(ambient::numeric::matrix<T>& a){
            a.fill_value(1.0);
        }
        
} } }

#undef ATOMIC
#endif
