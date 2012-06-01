#ifndef __ALPS_RESIZABLE_MATRIX_CONCEPT_CHECK_HPP__
#define __ALPS_RESIZABLE_MATRIX_CONCEPT_CHECK_HPP__
#include <boost/concept_check.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <stdexcept>
#include "types/dense_matrix/matrix_concept_check.hpp"

namespace maquis {
    namespace types {

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
            append_rows(x, std::make_pair(dataA.begin(),dataA.end()) );
            append_rows(x, std::make_pair(dataA.begin(),dataA.end()),1);
            append_cols(x, std::make_pair(dataB.begin(),dataB.end()) );
            append_cols(x, std::make_pair(dataB.begin(),dataB.end()),1);
            
            // Remove
            remove_rows(x,1);
            remove_rows(x,1,1);
            remove_cols(x,1);
            remove_cols(x,1,1);
    
            // Insert
            insert_rows(x,1, std::make_pair(dataA.begin(),dataA.end()) );
            insert_rows(x,1, std::make_pair(dataA.begin(),dataA.end()),1);
            insert_cols(x,1, std::make_pair(dataB.begin(),dataB.end()) );
            insert_cols(x,1, std::make_pair(dataB.begin(),dataB.end()),1); 
        }
    };

    }  // namespace type   
} // namespace maquis 
#endif //__ALPS_RESIZABLE_MATRIX_CONCEPT_CHECK_HPP__
