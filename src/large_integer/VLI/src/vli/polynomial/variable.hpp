#ifndef VLI_POLYNOMIAL_VARIABLE_HPP
#define VLI_POLYNOMIAL_VARIABLE_HPP
#include <boost/mpl/char.hpp>

namespace vli {

/*! \class var
        \brief This class models a variable

        This class just encapsulates the char of the variables of the polynomials
*/
template <char X>
class var
: public boost::mpl::char_<X> {
};

/* \cond I do not need this part in the doc*/
class no_variable {
};

template <int Order>
struct max_order_each {
    static unsigned int const value = Order;
};

template <int Order>
struct max_order_combined {
    static unsigned int const value = Order;
};

/* \endcond I do not need this part in the doc*/

} // end namespace vli
#endif //VLI_POLYNOMIAL_VARIABLE_HPP
