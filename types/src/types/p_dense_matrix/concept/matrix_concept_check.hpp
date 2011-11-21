#ifndef __ALPS_MATRIX_CONCEPT_CHECK_HPP__
#define __ALPS_MATRIX_CONCEPT_CHECK_HPP__
#include <boost/concept_check.hpp>
#include <boost/type_traits/remove_const.hpp>

namespace maquis {
    namespace types {
    
    template <typename X>
    struct Matrix
            : boost::Assignable<X>
    {
        public:
            typedef typename X::value_type                          value_type;
            // TODO BOOST_CONCEPT_ASSERT((Field<value_type>));
            typedef typename X::size_type                           size_type;
            BOOST_CONCEPT_ASSERT((boost::UnsignedInteger<size_type>));
            typedef typename X::difference_type                     difference_type;
            BOOST_CONCEPT_ASSERT((boost::SignedInteger<difference_type>));
            
        BOOST_CONCEPT_USAGE(Matrix)
        {
            // Constructor
            typename boost::remove_const<X>::type x(1,1);
            // Copy constructor
            const X y(x);
            typename boost::remove_const<X>::type z = x;
    
            // Swap
            std::swap(x,z);
    
            // num_rows(), num_cols()
            std::size_t s = num_rows(y);
            s = num_cols(y);
            
            // Element access
            t = x(0,0);
            x(0,0)+=value_type();
    
            // operators
            z = x;
            x += y;
            x -= y;
            x *= t;
    
            z = x + y;
            z = x - y;
            z = x * y;
    
            // Matrix vector multiplication
            // this does not check for mixed types
    //#warning FIXME
            /* Which vector class is this supposed to use in the general case? */
    //        vector<value_type> v;
    //        v = x * v;
    
        }
    
        private:
            // Default constructable value_type
            value_type t;
    };
    
    } //namespace types
}// namespace maquis    
#endif //__ALPS_MATRIX_CONCEPT_CHECK_HPP__
