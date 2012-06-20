#ifndef MATRIX_BINDINGS_H
#define MATRIX_BINDINGS_H

#include "alps/numeric/matrix/matrix.hpp"
#include "alps/numeric/matrix/matrix_interface.hpp"
#include "alps/numeric/diagonal_matrix.hpp"

namespace maquis{ namespace traits{

    template <typename O, typename I, typename T = typename I::value_type> 
    struct binding {
        typedef T value_type;
        static O convert(const I& m){
            return static_cast<O>(m); // default
        }
    };

    template <typename T>
    struct binding< std::vector<T>, alps::numeric::diagonal_matrix<T> > {
        static std::vector<T> convert(const alps::numeric::diagonal_matrix<T>& m){
            return m.get_values();
        }
    };

    template<typename O, typename I> 
    O matrix_cast(I const& input){
       return binding<O,I>::convert(input);
    }

} }

#ifdef AMBIENT 
#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/p_diagonal_matrix.h"

namespace maquis{ namespace traits{

    template <typename T>
    struct binding< maquis::types::p_dense_matrix<T>, alps::numeric::matrix<T> > {
        static maquis::types::p_dense_matrix<T> convert(const alps::numeric::matrix<T>& m){
            size_t num_rows = m.num_rows();
            size_t num_cols = m.num_cols();
            size_t lda = m.stride2();
            maquis::types::p_dense_matrix<T> pm(num_rows, num_cols);    
            const std::vector<typename alps::numeric::matrix<T>::value_type>* v_ptr = &m.get_values();
            ambient::push< ambient::cast_to_p_dense<T> >(v_ptr, pm, num_rows, num_cols, lda);
            ambient::playout();
            return pm;
        }
    };

    template <typename T>
    struct binding< alps::numeric::matrix<T>, maquis::types::p_dense_matrix<T> > {
        static alps::numeric::matrix<T> convert(const maquis::types::p_dense_matrix<T>& pm){
            size_t num_rows = pm.num_rows();
            size_t num_cols = pm.num_cols();
            alps::numeric::matrix<T> m(num_rows, num_cols);    
            std::vector<typename alps::numeric::matrix<T>::value_type>* v_ptr = &m.get_values();
            ambient::push< ambient::cast_to_dense<T> >(v_ptr, pm, num_rows, num_cols);
            ambient::playout();
            return m;
        }
    };

    template <typename T>
    struct binding< alps::numeric::matrix<T>, maquis::types::p_dense_matrix_impl<T> > {
        static alps::numeric::matrix<T> convert(const maquis::types::p_dense_matrix_impl<T>& pm){
            size_t num_rows = pm.num_rows();
            size_t num_cols = pm.num_cols();
            alps::numeric::matrix<T> m(num_rows, num_cols);    
            std::vector<typename alps::numeric::matrix<T>::value_type>* v_ptr = &m.get_values();
            ambient::push< ambient::cast_to_dense<T> >(v_ptr, pm, num_rows, num_cols);
            ambient::playout();
            return m;
        }
    };

    template <typename T>
    struct binding< maquis::types::p_diagonal_matrix<T>, alps::numeric::diagonal_matrix<T> > {
        static maquis::types::p_diagonal_matrix<T> convert(const alps::numeric::diagonal_matrix<T>& m){
            size_t num_rows(m.num_rows());
            size_t num_cols(1);
            maquis::types::p_diagonal_matrix<T> pm(num_rows, num_rows);    
            const std::vector<typename alps::numeric::diagonal_matrix<T>::value_type>* v_ptr = &m.get_values();
            ambient::push< ambient::cast_to_p_dense<T> >(v_ptr, pm, num_rows, num_cols, num_rows);
            ambient::playout();
            return pm;
        }
    };

    template <typename T>
    struct binding< alps::numeric::diagonal_matrix<T>, maquis::types::p_diagonal_matrix<T> > {
        static alps::numeric::diagonal_matrix<T> convert(const maquis::types::p_diagonal_matrix<T>& pm){
            size_t num_cols(1);
            size_t num_rows = pm.num_rows();
            alps::numeric::diagonal_matrix<T> m(num_rows, num_rows);    
            std::vector<typename alps::numeric::diagonal_matrix<T>::value_type>* v_ptr = &m.get_values();
            ambient::push< ambient::cast_to_dense<T> >(v_ptr, pm, num_rows, num_cols);
            ambient::playout();
            return m;
        }
    };

    template <typename T>
    struct binding< std::vector<T>, maquis::types::p_diagonal_matrix<T> > {
        static std::vector<T> convert(const maquis::types::p_diagonal_matrix<T>& pm){
            return binding<alps::numeric::diagonal_matrix<T>, maquis::types::p_diagonal_matrix<T> >::convert(pm).get_values();
        }
    };

} }

template<typename T>
bool operator == (alps::numeric::matrix<T> const & a, maquis::types::p_dense_matrix_impl<T> const & b) 
{
    ambient::future<int> ret(1);
    maquis::types::p_dense_matrix<T> pa = maquis::traits::matrix_cast<maquis::types::p_dense_matrix<T> >(a);
    ambient::push< ambient::validation<T> >(pa, b, ret); 
    return ((int)ret > 0);
}

template<typename T>
bool operator == (alps::numeric::matrix<T> const & a, maquis::types::p_dense_matrix<T> const & b) 
{
    ambient::future<int> ret(1);
    maquis::types::p_dense_matrix<T> pa = maquis::traits::matrix_cast<maquis::types::p_dense_matrix<T> >(a);
    ambient::push< ambient::validation<T> >(pa, b, ret); 
    return ((int)ret > 0);
}

template<typename T>
bool operator == (alps::numeric::diagonal_matrix<T> const & a, maquis::types::p_diagonal_matrix<T> const & b) 
{
    ambient::future<int> ret(1);
    maquis::types::p_diagonal_matrix<T> pa = maquis::traits::matrix_cast<maquis::types::p_diagonal_matrix<T> >(a);
    ambient::push< ambient::validation<T> >(pa, b, ret); 
    return ((int)ret > 0);
}

template<typename T>
bool operator == (maquis::types::p_dense_matrix<T> const & a, maquis::types::p_dense_matrix<T> const & b) 
{
    ambient::future<int> ret(1);
    ambient::push< ambient::validation<T> >(a, b, ret); 
    return ((int)ret > 0);
}

template<typename T>
bool operator == (maquis::types::p_dense_matrix_impl<T> const & pm, alps::numeric::matrix<T> const & m){
    return (m == pm);
}

template<typename T>
bool operator == (maquis::types::p_dense_matrix<T> const & pm, alps::numeric::matrix<T> const & m){
    return (m == pm);
}

template<typename T>
bool operator == (maquis::types::p_diagonal_matrix<T> const & pm, alps::numeric::diagonal_matrix<T> const & m){
    return (m == pm);
}

#endif
#endif
