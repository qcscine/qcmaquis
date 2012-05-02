#ifndef __ALPS_TYPES_TRAITS_HPP__
#define __ALPS_TYPES_TRAITS_HPP__

#include "utils/types.h"
#include "types/utils/matrix_vector_traits.h"
#include <alps/numeric/real.hpp>

namespace maquis {
    namespace types {

    template <typename T> class p_dense_matrix;
    template <typename T> class p_diagonal_matrix;
    
    }
}

namespace maquis {
    namespace traits {

    template<class T> struct scalar_type { typedef typename T::value_type type; };
    template<class T> struct scalar_type <maquis::types::p_dense_matrix<T> > { typedef typename maquis::types::p_dense_matrix<T>::scalar_type type; };
    template<class T> struct scalar_type <maquis::types::p_diagonal_matrix<T> > { typedef typename maquis::types::p_dense_matrix<T>::scalar_type type; };

    template <class T>
    inline typename utils::real_type<T>::type real(T x){
        return alps::numeric::real(x);
    }

    template <class T>
    inline std::vector<typename utils::real_type<T>::type> real(std::vector<T> x) 
    {
        std::vector<typename utils::real_type<T>::type> re; re.reserve(x.size());
        std::transform(x.begin(),x.end(),std::back_inserter(re),
                       static_cast<typename utils::real_type<T>::type (*)(T)>(&real));
        return re;
    }

    }
}

#endif
