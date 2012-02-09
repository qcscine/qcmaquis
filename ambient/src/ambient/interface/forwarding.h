#ifndef __AMBIENT_TYPES_FORWARDING_HPP__
#define __AMBIENT_TYPES_FORWARDING_HPP__

namespace maquis { namespace types {

    template <class T>
    class p_dense_matrix_impl;

} }

namespace ambient { namespace models {

    template <typename S> 
    struct info < maquis::types::p_dense_matrix_impl<S> > {
        typedef parallel_t_info< maquis::types::p_dense_matrix_impl<S> > typed; 
        typedef S value_type;
    };

    template <typename S> 
    struct info < const maquis::types::p_dense_matrix_impl<S> > { 
        typedef parallel_t_info< maquis::types::p_dense_matrix_impl<S> > typed; 
        typedef S value_type;
    };

} }

#endif
