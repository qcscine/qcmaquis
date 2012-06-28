#ifndef __ALPS_TYPES_TRAITS_HPP__
#define __ALPS_TYPES_TRAITS_HPP__

#include "utils/types.h"
#include <alps/numeric/real.hpp>
#include <alps/numeric/matrix/matrix_traits.hpp>

namespace maquis { namespace types {

    struct Transpose   { 
        static const char* code(){ return "T"; } 
        template<class Matrix> static typename Matrix::size_type first(Matrix const & m){ return m.num_cols(); }
        template<class Matrix> static typename Matrix::size_type second(Matrix const & m){ return m.num_rows(); }
        template<class Matrix> static Matrix eval(const Matrix& m){ return transpose(m); }
    };
    struct NoTranspose { 
        static const char* code(){ return "N"; } 
        template<class Matrix> static typename Matrix::size_type first(Matrix const & m){ return m.num_rows(); }
        template<class Matrix> static typename Matrix::size_type second(Matrix const & m){ return m.num_cols(); }
        template<class Matrix> static const Matrix& eval(const Matrix& m){ return m; }
    };
} }

namespace maquis { namespace traits {

    template<class T> struct scalar_type { typedef typename T::value_type type; };
   
} }

#ifdef USE_AMBIENT
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

    template<class T> struct scalar_type <ambient::numeric::matrix<T> > { typedef typename ambient::numeric::matrix<T>::scalar_type type; };
    template<class T> struct scalar_type <ambient::numeric::diagonal_matrix<T> > { typedef typename ambient::numeric::matrix<T>::scalar_type type; };

    template <class T> 
    inline typename utils::real_type<T>::type real(T x){
        return alps::numeric::real((typename T::value_type)x);
    }

    template <>
    inline utils::real_type<double>::type real(double x){
        return x;
    }

    template <typename T>
    inline typename utils::real_type<double>::type real(std::complex<T> x){
        return alps::numeric::real(x);
    }

    template <class T>
    inline std::vector<typename utils::real_type<T>::type> real(std::vector<T> x){
        std::vector<typename utils::real_type<T>::type> re; re.reserve(x.size());
        std::transform(x.begin(),x.end(),std::back_inserter(re),
                       static_cast<typename utils::real_type<T>::type (*)(T)>(&real));
        return re;
    }

} }
#endif

#endif
