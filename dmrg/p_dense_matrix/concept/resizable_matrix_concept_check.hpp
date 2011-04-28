#ifndef __ALPS_RESIZABLE_MATRIX_CONCEPT_CHECK_HPP__
#define __ALPS_RESIZABLE_MATRIX_CONCEPT_CHECK_HPP__
#include <boost/concept_check.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <stdexcept>
#include "p_dense_matrix/concept/matrix_concept_check.hpp"

namespace blas
{

template <typename X>
struct ResizableMatrix
        : Matrix<X>
{
    public:
    BOOST_CONCEPT_USAGE(ResizableMatrix)
    {
        typename boost::remove_const<X>::type x(1,1);

        // Resize
        resize(x,2,2);

        // Append
        std::vector<typename X::value_type> dataA(2,2);
        std::vector<typename X::value_type> dataB(4,2);
        
        // Remove
        remove_rows(x,1);
        remove_rows(x,1,1);
        remove_cols(x,1);
        remove_cols(x,1,1);
    }
};

}

#endif //__ALPS_RESIZABLE_MATRIX_CONCEPT_CHECK_HPP__
