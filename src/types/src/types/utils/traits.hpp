#ifndef __ALPS_TYPES_TRAITS_HPP__
#define __ALPS_TYPES_TRAITS_HPP__

#include "utils/types.h"
#include <alps/numeric/real.hpp>
#include <alps/numeric/matrix/matrix_traits.hpp>

namespace maquis { namespace types {

    template<typename T> class p_dense_matrix;
    template<typename T> class p_diagonal_matrix;

    class Transpose { };
    class NoTranspose { };
    
    namespace detail {

        template<class Tag> struct evaluate_tag;
        template<> struct evaluate_tag<maquis::types::NoTranspose> {
            template<class Matrix>
            static Matrix const & eval(Matrix const & m) { return m; }
        };
        template<> struct evaluate_tag<maquis::types::Transpose> {
            template<class Matrix>
            static Matrix eval(Matrix const & m) { return transpose(m); }
        };
        
        template<class Tag> struct dims;
        template<> struct dims<NoTranspose> {
            template<class Matrix> 
            static typename Matrix::size_type first(Matrix const & m){ return num_rows(m); }
            template<class Matrix>
            static typename Matrix::size_type second(Matrix const & m){ return num_cols(m); }
        };
        template<> struct dims<Transpose> {
            template<class Matrix>
            static typename Matrix::size_type first(Matrix const & m){ return num_cols(m); }
            template<class Matrix>
            static typename Matrix::size_type second(Matrix const & m){ return num_rows(m); }
        };
    }

} } // namespace maquis::types

namespace alps { namespace numeric {

    template<typename T> 
    struct associated_diagonal_matrix< maquis::types::p_dense_matrix<T> > {
        typedef maquis::types::p_diagonal_matrix<T> type;
    };
    template<typename T> 
    struct associated_real_diagonal_matrix< maquis::types::p_dense_matrix<T> > {
        typedef maquis::types::p_diagonal_matrix<typename utils::real_type<T>::type> type;
    };
    template<typename T> 
    struct associated_vector< maquis::types::p_dense_matrix<T> > {
        typedef std::vector<T> type;
    };
    template<typename T> 
    struct associated_real_vector< maquis::types::p_dense_matrix<T> > {
        typedef std::vector<typename utils::real_type<T>::type> type;
    };

} }

namespace maquis { namespace traits {

    template<class T> struct scalar_type { typedef typename T::value_type type; };
    template<class T> struct scalar_type <maquis::types::p_dense_matrix<T> > { typedef typename maquis::types::p_dense_matrix<T>::scalar_type type; };
    template<class T> struct scalar_type <maquis::types::p_diagonal_matrix<T> > { typedef typename maquis::types::p_dense_matrix<T>::scalar_type type; };

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

template<class Matrix1, class Matrix2, class Matrix3, class Tag1, class Tag2>
void gemm(Matrix1 const & A, Tag1, Matrix2 const & B, Tag2, Matrix3 & C){
    gemm(maquis::types::detail::evaluate_tag<Tag1>::eval(A),
         maquis::types::detail::evaluate_tag<Tag2>::eval(B),
         C);
}

template<class Matrix1, class Matrix2, class Tag1, class Tag2>
std::pair<std::size_t, std::size_t> result_size(Matrix1 const & A, Tag1, Matrix2 const & B, Tag2){
    return std::make_pair(maquis::types::detail::dims<Tag1>::first(A), maquis::types::detail::dims<Tag2>::second(B));
}

#endif
