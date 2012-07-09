#ifndef __ALPS_TYPES_TRAITS_HPP__
#define __ALPS_TYPES_TRAITS_HPP__

#include "utils/types.h"
#include <alps/numeric/real.hpp>
#include <alps/numeric/matrix/matrix_traits.hpp>
#include <boost/numeric/bindings/tag.hpp>

namespace ambient { namespace numeric {

    template<typename T> class matrix;
    template<typename T> class diagonal_matrix;

} }

namespace maquis { namespace types {

    struct Transpose   { 
        typedef boost::numeric::bindings::tag::transpose type;
        template<class Matrix> static typename Matrix::size_type first(Matrix const & m){ return m.num_cols(); }
        template<class Matrix> static typename Matrix::size_type second(Matrix const & m){ return m.num_rows(); }
        template<class Matrix> static Matrix eval(const Matrix& m){ return transpose(m); }
        static const char* code(){ return "T"; } 
    };
    struct NoTranspose { 
        typedef boost::numeric::bindings::tag::no_transpose type;
        template<class Matrix> static typename Matrix::size_type first(Matrix const & m){ return m.num_rows(); }
        template<class Matrix> static typename Matrix::size_type second(Matrix const & m){ return m.num_cols(); }
        template<class Matrix> static const Matrix& eval(const Matrix& m){ return m; }
        static const char* code(){ return "N"; } 
    };

} }

namespace alps { namespace numeric {

    template<typename T> 
    struct associated_diagonal_matrix< ambient::numeric::matrix<T> > {
        typedef ambient::numeric::diagonal_matrix<T> type;
    };
    template<typename T> 
    struct associated_real_diagonal_matrix< ambient::numeric::matrix<T> > {
        typedef ambient::numeric::diagonal_matrix<typename utils::real_type<T>::type> type;
    };
    template<typename T> 
    struct associated_vector< ambient::numeric::matrix<T> > {
        typedef std::vector<T> type;
    };
    template<typename T> 
    struct associated_real_vector< ambient::numeric::matrix<T> > {
        typedef std::vector<typename utils::real_type<T>::type> type;
    };

} }

namespace maquis { namespace traits {

    template<class T> struct scalar_type { typedef typename T::value_type type; };
    template<class T> struct scalar_type <ambient::numeric::matrix<T> > { typedef typename ambient::numeric::matrix<T>::scalar_type type; };
    template<class T> struct scalar_type <ambient::numeric::diagonal_matrix<T> > { typedef typename ambient::numeric::matrix<T>::scalar_type type; };

} }

#endif
