/* \cond I do not need this part in the doc*/

#ifndef VLI_POLYNOMIAL_DETAIL_HELPERS_HPP
#define VLI_POLYNOMIAL_DETAIL_HELPERS_HPP

namespace vli {
namespace detail {
namespace max_order_combined_helpers {

    // TODO reuse a different factorial...
    template <unsigned int N>
    struct factorial {
        static unsigned int const value = N*factorial<N-1>::value;
    };

    template <>
    struct factorial<0> {
        static unsigned int const value = 1;
    };


    template<unsigned int NK, unsigned int K>
    struct size_helper {
        static unsigned int const value = NK*size_helper<NK-1,K>::value;
    };

    template <unsigned int K>
    struct size_helper<K,K> {
        static unsigned int const value = 1;
    };

    template <unsigned int N, unsigned int K>
    struct size {
        // N variables, max order K -> n+k-1 over k  = (n+k-1)! / ( (n-1)! k! ) combinations
        // Assuming N > 0
        static unsigned int const value = size_helper<N+K-1,K>::value/factorial<N-1>::value;
    };

}

template <int Var, int NumVars, int Order>
struct stride {
    //Var+1???
    static unsigned int const value = Var < NumVars ? Order+1 : 1;
};

template <class Var0, class Var1, class Var2, class Var3> 	 	 
struct num_of_variables_helper { 	 	 
    // May be generalized using boost MPL 	 	 
    static unsigned int const value = 4; 	 	 
}; 	 	 
 	 	 
template <class Var0, class Var1, class Var2> 	 	 
struct num_of_variables_helper<Var0, Var1, Var2, no_variable> { 	 	 
    static unsigned int const value = 3; 	 	 
}; 	 	 
     	 
template <class Var0, class Var1> 	 	 
struct num_of_variables_helper<Var0, Var1, no_variable, no_variable> { 	 	 
    static unsigned int const value = 2; 	 	 
}; 	 	 
     	 
template <class Var0> 	 	 
struct num_of_variables_helper<Var0, no_variable, no_variable, no_variable> { 	 	 
    static unsigned int const value = 1; 	 	 
}; 	 	 

template <class MaxOrder, int NumVars, int Coeff=1> // this coeff is equal to 2 when I start the thread grid 
struct num_coefficients;

template <int Order, int NumVars, int Coeff>
struct num_coefficients<max_order_each<Order>, NumVars, Coeff>{
    static unsigned int const value = (Coeff*Order+1)*num_coefficients<max_order_each<Order>, NumVars-1, Coeff>::value;
};

template <int Order, int Coeff>
struct num_coefficients<max_order_each<Order>, 0, Coeff>{
    static unsigned int const value = 1;
};

template<int Order, int NumVars, int Coeff>
struct num_coefficients<max_order_combined<Order>, NumVars, Coeff>{
    static unsigned int const value = detail::max_order_combined_helpers::size<NumVars+1, Coeff*Order>::value;
};

} // end namespace detail
} // end namespace vli
#endif //VLI_POLYNOMIAL_DETAIL_HELPERS_HPP

/* \endcond I do not need this part in the doc*/