#ifndef VLI_POLYNOMIAL_VARIABLE_HPP
#define VLI_POLYNOMIAL_VARIABLE_HPP
#include <boost/mpl/char.hpp>

namespace vli {

template <char X>
class var
: public boost::mpl::char_<X> {
};

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

} // end namespace vli
#endif //VLI_POLYNOMIAL_VARIABLE_HPP
