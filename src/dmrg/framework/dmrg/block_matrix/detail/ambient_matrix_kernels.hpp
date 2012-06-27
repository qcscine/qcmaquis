/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2012 by Alexandr Kosenkov <alex.kosenkov@gmail.com>
 *                            Timothee Ewart <timothee.ewart@gmail.com>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_KERNELS_P_DENSE_MATRIX_HPP
#define MAQUIS_DMRG_KERNELS_P_DENSE_MATRIX_HPP

#include <types/p_dense_matrix/p_dense_matrix.h>

namespace maquis { namespace dmrg { namespace detail {
        
        template <typename T>
        inline void reshape_l2r(const maquis::types::p_dense_matrix<T>& left, maquis::types::p_dense_matrix<T>& right,
                                size_t left_offset, size_t right_offset, 
                                size_t sdim, size_t ldim, size_t rdim)
        { // gs
#ifdef AMBIENT_SERIAL_CHECK
            alps::numeric::matrix<T> sl = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(left);
            alps::numeric::matrix<T> sr = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(right);
            reshape_l2r(sl, sr, left_offset, right_offset, sdim, ldim, rdim);
#endif
            USE_ATOMIC(left.atomic() && right.atomic(), reshape_l2r, 
                       left, right, left_offset, right_offset, sdim, ldim, rdim);
#ifdef AMBIENT_SERIAL_CHECK
            if(sr == right){} else printf("--------------------- RESHAPE L2R WAS INCORRECT!\n");
#endif
        }
        
        template <typename T>
        inline void reshape_r2l(maquis::types::p_dense_matrix<T>& left, const maquis::types::p_dense_matrix<T>& right,
                                size_t left_offset, size_t right_offset, 
                                size_t sdim, size_t ldim, size_t rdim)
        { // gs
#ifdef AMBIENT_SERIAL_CHECK
            alps::numeric::matrix<T> sl = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(left);
            alps::numeric::matrix<T> sr = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(right);
            reshape_r2l(sl, sr, left_offset, right_offset, sdim, ldim, rdim);
#endif
            USE_ATOMIC(left.atomic() && right.atomic(), reshape_r2l, 
                       left, right, left_offset, right_offset, sdim, ldim, rdim);
#ifdef AMBIENT_SERIAL_CHECK
            if(sl == left){} else printf("--------------------- RESHAPE R2L WAS INCORRECT!\n");
#endif
        }
        
        template <typename T>
        inline void lb_tensor_mpo(maquis::types::p_dense_matrix<T>& out, const maquis::types::p_dense_matrix<T>& in, const maquis::types::p_dense_matrix<T>& alfa,
                                  size_t out_offset, size_t in_offset, 
                                  size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
        { // gs
#ifdef AMBIENT_SERIAL_CHECK
            alps::numeric::matrix<T> sout = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(out);
            alps::numeric::matrix<T> sin = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(in);
            alps::numeric::matrix<T> salfa = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(alfa);
            lb_tensor_mpo(sout, sin, salfa, out_offset, in_offset, sdim1, sdim2, ldim, rdim);
#endif
            USE_ATOMIC(out.atomic() && in.atomic() && alfa.atomic(), lb_tensor_mpo, 
                       out, in, alfa, out_offset, in_offset, sdim1, sdim2, ldim, rdim);
#ifdef AMBIENT_SERIAL_CHECK
            if(sout == out){} else printf("--------------------- LB TENSOR MPO WAS INCORRECT!\n");
#endif
        }
        
        template <typename T>
        inline void rb_tensor_mpo(maquis::types::p_dense_matrix<T>& out, const maquis::types::p_dense_matrix<T>& in, const maquis::types::p_dense_matrix<T>& alfa,
                                  size_t out_offset, size_t in_offset, 
                                  size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
        { // gs
#ifdef AMBIENT_SERIAL_CHECK
            alps::numeric::matrix<T> sout = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(out);
            alps::numeric::matrix<T> sin = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(in);
            alps::numeric::matrix<T> salfa = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(alfa);
            rb_tensor_mpo(sout, sin, salfa, out_offset, in_offset, sdim1, sdim2, ldim, rdim);
#endif
            USE_ATOMIC(out.atomic() && in.atomic() && alfa.atomic(), rb_tensor_mpo, 
                       out, in, alfa, out_offset, in_offset, sdim1, sdim2, ldim, rdim);
#ifdef AMBIENT_SERIAL_CHECK
            if(sout == out){} else printf("--------------------- RB TENSOR MPO WAS INCORRECT!\n");
#endif
        }
        
        template <typename T>
        inline void norm(const maquis::types::p_dense_matrix<T>& a, scalar_type& ret){
            // gs
#ifdef AMBIENT_SERIAL_CHECK
            alps::numeric::matrix<T> s = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(a);
            typename alps::numeric::matrix<T>::value_type sret((T)ret);
            norm(s, sret);
#endif
            USE_ATOMIC(a.atomic(), scalar_norm, a, ret);
#ifdef AMBIENT_SERIAL_CHECK
            if(std::abs((T)ret - sret) > 0.01) printf("--------------------- SCALAR NORM IS INCORRECT (%.2f vs %.2f)\n", (T)ret, sret); 
#endif
        }
        
        template <typename T>
        inline void overlap(const maquis::types::p_dense_matrix<T>& a, const maquis::types::p_dense_matrix<T>& b, scalar_type& ret){
            // gs
#ifdef AMBIENT_SERIAL_CHECK
            alps::numeric::matrix<T> sa = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(a);
            alps::numeric::matrix<T> sb = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(b);
            T sret((T)ret); overlap(sa, sb, sret);
#endif
            USE_ATOMIC(a.atomic(), overlap, a, b, ret);
#ifdef AMBIENT_SERIAL_CHECK
            if(std::abs((T)ret - sret) > 0.01) printf("--------------------- SCALAR NORM OVERLAP IS INCORRECT (%.2f vs %.2f)\n", (T)ret, sret); 
#endif
        }
        
        template <typename T>
        inline void bond_renyi_entropies(const maquis::types::p_diagonal_matrix<T>& m, typename alps::numeric::associated_real_vector<maquis::types::p_dense_matrix<T> >::type& sv){
            // gs
#ifdef AMBIENT_SERIAL_CHECK
            alps::numeric::diagonal_matrix<T> sm = maquis::traits::matrix_cast<alps::numeric::diagonal_matrix<T> >(m);
            typename alps::numeric::associated_real_vector<maquis::types::p_dense_matrix<T> >::type ssv(sv);
            bond_renyi_entropies(sm, ssv);
#endif
            std::vector<T>* sc_ptr = &sv;
            ambient::push< ambient::push_back_sqr_gt<T> >(m, sc_ptr);
            ambient::playout();
#ifdef AMBIENT_SERIAL_CHECK
            if(sv != ssv) printf("--------------------- BOND RENYI ENTROPIES WAS INCORRECT!\n");
#endif
        }
        
        template <typename T>
        inline void left_right_boundary_init(maquis::types::p_dense_matrix<T>& a){
            // gs
#ifdef AMBIENT_SERIAL_CHECK
            alps::numeric::matrix<T> sa = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(a);
            left_right_boundary_init(sa);
#endif
            a.fill_value(1.0);
#ifdef AMBIENT_SERIAL_CHECK
            if(sa == a){}else printf("--------------------- INCORRECT LEFT RIGHT BOUNDARY INIT!\n");
#endif
        }
        
} } } // namespace maquis::dmrg::detail


#endif // MAQUIS_DMRG_KERNELS_P_DENSE_MATRIX_HPP
