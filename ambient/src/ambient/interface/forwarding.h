#ifndef AMBIENT_TYPES_FORWARDING
#define AMBIENT_TYPES_FORWARDING

namespace maquis { namespace types {

    template <class T>
    class p_dense_matrix;

    template <class T>
    class p_diagonal_matrix;

    template <class T>
    class p_dense_matrix_impl;

} }

namespace ambient { 

    template <typename T>
    class future;

    #define BOOST_SP_NO_SP_CONVERTIBLE
    #include <boost/intrusive_ptr.hpp>
 
    template <typename T>
    void copy_l(T& ac, pinned const T& a);

    template <typename T>
    void copy_c(T& ac, pinned const T& a);

namespace models {

    template <typename T> 
    struct info { 
        typedef singular_t_info<T> typed; 
        static inline T& unfold(T& naked){ return naked; }
    };

    template <typename T>
    struct info < boost::intrusive_ptr<T> > { 
    };

    template <typename T>
    struct info < const boost::intrusive_ptr<T> > { 
    };

    template <typename S>
    struct info < ambient::future<S> > { 
        typedef ambient::future<S> T;
        static inline S*& unfold(T& folded){
            return folded.unfold();
        }
    };

    template <typename S>
    struct info < const ambient::future<S> > { 
        typedef const ambient::future<S> T;
        static inline const S*& unfold(T& folded){
            assert(false); // remove if believe in const future
            return folded.unfold();
        }
    };

    template <typename S>
    struct info < maquis::types::p_diagonal_matrix<S> > { 
        typedef maquis::types::p_diagonal_matrix<S> T;
        static inline maquis::types::p_dense_matrix_impl<S>& unfold(T& folded){
            return *folded.get_data().impl;
        }
    };

    template <typename S>
    struct info < const maquis::types::p_diagonal_matrix<S> > { 
        typedef const maquis::types::p_diagonal_matrix<S> T;
        static inline const maquis::types::p_dense_matrix_impl<S>& unfold(T& folded){
            return *folded.get_data().impl;
        }
    };

    template <typename S>
    struct info < maquis::types::p_dense_matrix<S> > { 
        typedef maquis::types::p_dense_matrix<S> T;
        static inline maquis::types::p_dense_matrix_impl<S>& unfold(T& folded){
            return *folded.impl;
        }
    };

    template <typename S>
    struct info < const maquis::types::p_dense_matrix<S> > { 
        typedef const maquis::types::p_dense_matrix<S> T;
        static inline const maquis::types::p_dense_matrix_impl<S>& unfold(T& folded){
            return *folded.impl;
        }
    };

    template <typename S>
    struct info < maquis::types::p_dense_matrix_impl<S> > {
        typedef maquis::types::p_dense_matrix_impl<S> T;
        static inline T& unfold(T& naked){ return naked; }
        typedef parallel_t_info< T > typed; 
        typedef S value_type;
    };

    template <typename S>
    struct info < const maquis::types::p_dense_matrix_impl<S> > { 
        typedef maquis::types::p_dense_matrix_impl<S> T;
        static inline const T& unfold(const T& naked){ return naked; }
        typedef parallel_t_info< T > typed; 
        typedef S value_type;
    };

} }

#endif
