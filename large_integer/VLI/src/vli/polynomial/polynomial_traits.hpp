#ifndef VLI__POLYNOMIAL_TRAITS_HPP
#define VLI__POLYNOMIAL_TRAITS_HPP
namespace vli {


class no_variable;

template <class Coeff, class OrderSpecification, class Var0, class Var1, class Var2, class Var3>
class polynomial;



template <typename Polynomial>
struct polynomial_multiply_result_type {
};

template <typename Coeff, unsigned int Order, class Var0, class Var1, class Var2, class Var3>
struct polynomial_multiply_result_type<polynomial<Coeff,max_order_each<Order>,Var0,Var1,Var2,Var3> > {
    typedef polynomial<Coeff,max_order_each<2*Order>,Var0,Var1,Var2,Var3> type;
};

template <typename Coeff, unsigned int Order, class Var0, class Var1, class Var2, class Var3>
struct polynomial_multiply_result_type<polynomial<Coeff,max_order_combined<Order>,Var0,Var1,Var2,Var3> > {
    typedef polynomial<Coeff,max_order_combined<2*Order>,Var0,Var1,Var2,Var3> type;
};



template <class Polynomial>
struct polynomial_multiply_keep_order_result_type {
};

template <typename Coeff, class OrderSpecification, class Var0, class Var1, class Var2, class Var3>
struct polynomial_multiply_keep_order_result_type< polynomial<Coeff,OrderSpecification,Var0,Var1,Var2,Var3> > {
    typedef polynomial<Coeff,OrderSpecification,Var0,Var1,Var2,Var3> type;
};


template <typename Polynomial>
struct exponent_type{
    typedef typename Polynomial::exponent_type type;
};

template <class Polynomial>
struct iterator {
    typedef typename Polynomial::iterator type;
};

template <class Polynomial>
struct const_iterator {
    typedef typename Polynomial::const_iterator type;
};

template <class Polynomial>
struct element_descriptor {
    typedef typename Polynomial::element_descriptor type;
};

template <class Polynomial>
struct order_specification {
    typedef typename Polynomial::order_specification type;
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

#endif //VLI__POLYNOMIAL_TRAITS_HPP
