#ifndef __AMBIENT_TYPES_FORWARDING_HPP__
#define __AMBIENT_TYPES_FORWARDING_HPP__

namespace maquis { namespace types {

    template <class T>
    class p_dense_matrix_impl;

    template <class T>
    class p_dense_matrix;

} }

namespace ambient { 

    #define BOOST_SP_NO_SP_CONVERTIBLE
    #include <boost/intrusive_ptr.hpp>
    
namespace models {

    template <typename T> 
    struct info { 
        typedef singular_t_info<T> typed; 
    };

    template <typename T>
    struct info < boost::intrusive_ptr<T> > { };

    template <typename T>
    struct info < const boost::intrusive_ptr<T> > { };

    template <typename S>
    struct info < maquis::types::p_dense_matrix<S> > { };

    template <typename S>
    struct info < const maquis::types::p_dense_matrix<S> > { };

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
