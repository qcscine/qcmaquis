#ifndef VLI__POLYNOMIAL_TRAITS_HPP
#define VLI__POLYNOMIAL_TRAITS_HPP

#include <boost/mpl/int.hpp>

/* \cond I do not need this part in the doc*/


namespace vli {

class no_variable;

template <class Coeff, class MaxOrder, class Var0, class Var1, class Var2, class Var3>
class polynomial;

template <typename Polynomial>
struct polynomial_multiply_result_type {
};

template <typename Coeff, int Order, class Var0, class Var1, class Var2, class Var3>
struct polynomial_multiply_result_type<polynomial<Coeff,max_order_each<Order>,Var0,Var1,Var2,Var3> > {
    typedef polynomial<Coeff,max_order_each<2*Order>,Var0,Var1,Var2,Var3> type;
};

template <typename Coeff, int Order, class Var0, class Var1, class Var2, class Var3>
struct polynomial_multiply_result_type<polynomial<Coeff,max_order_combined<Order>,Var0,Var1,Var2,Var3> > {
    typedef polynomial<Coeff,max_order_combined<2*Order>,Var0,Var1,Var2,Var3> type;
};

template <class Polynomial>
struct polynomial_multiply_keep_order_result_type {
};

template <typename Coeff, class MaxOrder, class Var0, class Var1, class Var2, class Var3>
struct polynomial_multiply_keep_order_result_type< polynomial<Coeff,MaxOrder,Var0,Var1,Var2,Var3> > {
    typedef polynomial<Coeff,MaxOrder,Var0,Var1,Var2,Var3> type;
};


template <class Polynomial>
struct max_order {
    typedef typename Polynomial::max_order type;
};

template <class Polynomial>
struct num_variables {
};

template <class Coeff, class MaxOrder, class Var0, class Var1, class Var2, class Var3>
struct num_variables<polynomial<Coeff,MaxOrder,Var0,Var1,Var2,Var3> >
: boost::mpl::int_<4> {
};

template <class Coeff, class MaxOrder, class Var0, class Var1, class Var2>
struct num_variables<polynomial<Coeff,MaxOrder,Var0,Var1,Var2,no_variable> >
: boost::mpl::int_<3> {
};

template <class Coeff, class MaxOrder, class Var0, class Var1>
struct num_variables<polynomial<Coeff,MaxOrder,Var0,Var1,no_variable,no_variable> >
: boost::mpl::int_<2> {
};

template <class Coeff, class MaxOrder, class Var0>
struct num_variables<polynomial<Coeff,MaxOrder,Var0,no_variable,no_variable,no_variable> >
: boost::mpl::int_<1> {
};

template <class Polynomial, unsigned int N>
struct variable {
    typedef no_variable type;
};

template <class Polynomial>
struct variable<Polynomial,0> {
    typedef typename Polynomial::var0_type type;
};

template <class Polynomial>
struct variable<Polynomial,1> {
    typedef typename Polynomial::var1_type type;
};

template <class Polynomial>
struct variable<Polynomial,2> {
    typedef typename Polynomial::var2_type type;
};

template <class Polynomial>
struct variable<Polynomial,3> {
    typedef typename Polynomial::var3_type type;
};

} // end namespace vli

/* \endcond I do not need this part in the doc*/


#endif //VLI__POLYNOMIAL_TRAITS_HPP
