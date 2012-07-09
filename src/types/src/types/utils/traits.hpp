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

namespace maquis { namespace traits {

    template<class T> struct scalar_type { typedef typename T::value_type type; };
    template<class T> struct scalar_type <ambient::numeric::matrix<T> > { typedef typename ambient::numeric::matrix<T>::scalar_type type; };
    template<class T> struct scalar_type <ambient::numeric::diagonal_matrix<T> > { typedef typename ambient::numeric::matrix<T>::scalar_type type; };

    template<class T> struct real_type { typedef typename real_type<typename T::value_type>::type type; };
    template<>        struct real_type<double> { typedef double type; };
    template<class T> struct real_type<std::complex<T> > { typedef T type; };

    template<class T> struct real_identity { static const T value = 1; };
    template<class T> struct imag_identity { static const T value = 1; };
    template<class T> struct real_identity<std::complex<T> > { static const std::complex<T> value; };
    template<class T> struct imag_identity<std::complex<T> > { static const std::complex<T> value; };
    template<class T> const std::complex<T> real_identity<std::complex<T> >::value = std::complex<T>(1,0);
    template<class T> const std::complex<T> imag_identity<std::complex<T> >::value = std::complex<T>(0,1);

} }

namespace alps { namespace numeric {

    template<typename T> 
    struct associated_diagonal_matrix< ambient::numeric::matrix<T> > {
        typedef ambient::numeric::diagonal_matrix<T> type;
    };
    template<typename T> 
    struct associated_real_diagonal_matrix< ambient::numeric::matrix<T> > {
        typedef ambient::numeric::diagonal_matrix<typename maquis::traits::real_type<T>::type> type;
    };
    template<typename T> 
    struct associated_vector< ambient::numeric::matrix<T> > {
        typedef std::vector<T> type;
    };
    template<typename T> 
    struct associated_real_vector< ambient::numeric::matrix<T> > {
        typedef std::vector<typename maquis::traits::real_type<T>::type> type;
    };

} }

#endif
