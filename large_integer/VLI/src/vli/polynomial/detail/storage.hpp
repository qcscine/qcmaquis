#ifndef VLI_POLYNOMIAL_STORAGE_HPP
#define VLI_POLYNOMIAL_STORAGE_HPP

#include <boost/array.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/contains.hpp>
#include <cassert>
#include <boost/static_assert.hpp>

namespace vli {
namespace detail {

    template <class Var>
    struct var_exponent {
        typedef Var var_type;
        var_exponent(unsigned int exp)
        : exp(exp){
        }
        unsigned int exp;
    };

    template <>
    struct var_exponent<no_variable> {
        var_exponent(unsigned int i) {
        }
        static unsigned int const exp = 0;
    };

    template <class Var>
    inline std::ostream& operator << (std::ostream& os, var_exponent<Var> const& ve) {
        if(ve.exp != 0)
            os<<"*"<<var_exponent<Var>::var_type::value<<"^"<<ve.exp;
        return os;
    }

    template <>
    inline std::ostream& operator << (std::ostream& os, var_exponent<no_variable> const& ve) {
        return os;
    }

    // Some forward declaration
    template <class Var0, class Var1, class Var2, class Var3>
    struct element_descriptor_impl;
    
    template <class Var0, class Var1, class Var2, class Var3, class VarN>
    static inline unsigned int const exponent(element_descriptor_impl<Var0,Var1,Var2,Var3> const& e, VarN dummy) {
        return 0;
    }
    
    template <class Var0, class Var1, class Var2, class Var3>
    inline unsigned int & exponent(element_descriptor_impl<Var0,Var1,Var2,Var3> & e, Var0 dummy) {
        return e.var0.exp;
    }
    
    template <class Var0, class Var1, class Var2, class Var3>
    inline unsigned int const& exponent(element_descriptor_impl<Var0,Var1,Var2,Var3> const& e, Var0 dummy) {
        return e.var0.exp;
    }
    
    template <class Var0, class Var1, class Var2, class Var3>
    inline unsigned int & exponent(element_descriptor_impl<Var0,Var1,Var2,Var3> & e, typename boost::disable_if<boost::is_same<Var1,no_variable>,Var1>::type dummy) {
        return e.var1.exp;
    }
    
    template <class Var0, class Var1, class Var2, class Var3>
    inline unsigned int const& exponent(element_descriptor_impl<Var0,Var1,Var2,Var3> const& e, typename boost::disable_if<boost::is_same<Var1,no_variable>,Var1>::type dummy) {
        return e.var1.exp;
    }
    
    template <class Var0, class Var1, class Var2, class Var3>
    inline unsigned int & exponent(element_descriptor_impl<Var0,Var1,Var2,Var3> & e, typename boost::disable_if<boost::is_same<Var2,no_variable>,Var2>::type dummy) {
        return e.var2.exp;
    }
    
    template <class Var0, class Var1, class Var2, class Var3>
    inline unsigned int const& exponent(element_descriptor_impl<Var0,Var1,Var2,Var3> const& e, typename boost::disable_if<boost::is_same<Var2,no_variable>,Var2>::type dummy) {
        return e.var2.exp;
    }
    
    template <class Var0, class Var1, class Var2, class Var3>
    inline unsigned int & exponent(element_descriptor_impl<Var0,Var1,Var2,Var3> & e, typename boost::disable_if<boost::is_same<Var3,no_variable>,Var3>::type dummy) {
        return e.var3.exp;
    }

    template <class Var0, class Var1, class Var2, class Var3>
    inline unsigned int const& exponent(element_descriptor_impl<Var0,Var1,Var2,Var3> const& e, typename boost::disable_if<boost::is_same<Var3,no_variable>,Var3>::type dummy) {
        return e.var3.exp;
    }

    
    
    template <class Var0, class Var1, class Var2, class Var3>
    struct element_descriptor_impl {
        explicit element_descriptor_impl(unsigned int exp0, unsigned int exp1 = 0, unsigned int exp2 = 0, unsigned int exp3 = 0)
        : var0(exp0), var1(exp1), var2(exp2), var3(exp3) {
        }

        template <class MVar0, class MVar1, class MVar2, class MVar3>
        explicit element_descriptor_impl(element_descriptor_impl<MVar0,MVar1,MVar2,MVar3> const& e)
        : var0(exponent(e,Var0())), var1(exponent(e,Var1())), var2(exponent(e,Var2())), var3(exponent(e,Var3())) {
            // The old element descriptor variables must be a subset of the new ones 
            typedef boost::mpl::vector<Var0,Var1,Var2,Var3,no_variable> arg_variables;
            BOOST_STATIC_ASSERT(( boost::mpl::contains<arg_variables,MVar0>::value ));
            BOOST_STATIC_ASSERT(( boost::mpl::contains<arg_variables,MVar1>::value ));
            BOOST_STATIC_ASSERT(( boost::mpl::contains<arg_variables,MVar2>::value ));
            BOOST_STATIC_ASSERT(( boost::mpl::contains<arg_variables,MVar3>::value ));
        }
        var_exponent<Var0> var0;
        var_exponent<Var1> var1;
        var_exponent<Var2> var2;
        var_exponent<Var3> var3;
    };
    
    struct element_descr_minus_assign_impl {
        template <class Var0A, class Var1A, class Var2A, class Var3A, class Var0B, class Var1B, class Var2B, class Var3B, class VarN>
        static void apply(element_descriptor_impl<Var0A,Var1A,Var2A,Var3A>& a, element_descriptor_impl<Var0B,Var1B,Var2B,Var3B> const& b, VarN dummy) {
            assert(exponent(a,VarN()) >= exponent(b,VarN()));
            exponent(a,VarN()) -= exponent(b,VarN());
        }

        template <class Var0A, class Var1A, class Var2A, class Var3A, class Var0B, class Var1B, class Var2B, class Var3B>
        static void apply(element_descriptor_impl<Var0A,Var1A,Var2A,Var3A>& a, element_descriptor_impl<Var0B,Var1B,Var2B,Var3B> const& b, no_variable dummy) {
        }
    };

    struct element_descr_plus_assign_impl {
        template <class Var0A, class Var1A, class Var2A, class Var3A, class Var0B, class Var1B, class Var2B, class Var3B, class VarN>
        static void apply(element_descriptor_impl<Var0A,Var1A,Var2A,Var3A>& a, element_descriptor_impl<Var0B,Var1B,Var2B,Var3B> const& b, VarN dummy) {
            exponent(a,VarN()) += exponent(b,VarN());
        }

        template <class Var0A, class Var1A, class Var2A, class Var3A, class Var0B, class Var1B, class Var2B, class Var3B>
        static void apply(element_descriptor_impl<Var0A,Var1A,Var2A,Var3A>& a, element_descriptor_impl<Var0B,Var1B,Var2B,Var3B> const& b, no_variable dummy) {
        }
    };
    
    template <class Var0A, class Var1A, class Var2A, class Var3A, class Var0B, class Var1B, class Var2B, class Var3B>
    element_descriptor_impl<Var0A,Var1A,Var2A,Var3A> operator + (element_descriptor_impl<Var0A,Var1A,Var2A,Var3A> a, element_descriptor_impl<Var0B,Var1B,Var2B,Var3B> const& b) {
        // TODO static assert that lhs Vars are a superset of rhs Vars
        element_descr_plus_assign_impl::apply(a,b,Var0A());
        element_descr_plus_assign_impl::apply(a,b,Var1A());
        element_descr_plus_assign_impl::apply(a,b,Var2A());
        element_descr_plus_assign_impl::apply(a,b,Var3A());
        return a;
    }
    
    template <class Var0A, class Var1A, class Var2A, class Var3A, class Var0B, class Var1B, class Var2B, class Var3B>
    element_descriptor_impl<Var0A,Var1A,Var2A,Var3A> operator - (element_descriptor_impl<Var0A,Var1A,Var2A,Var3A> a, element_descriptor_impl<Var0B,Var1B,Var2B,Var3B> const& b) {
        // TODO static assert that lhs Vars are a superset of rhs Vars
        element_descr_minus_assign_impl::apply(a,b,Var0A());
        element_descr_minus_assign_impl::apply(a,b,Var1A());
        element_descr_minus_assign_impl::apply(a,b,Var2A());
        element_descr_minus_assign_impl::apply(a,b,Var3A());
        return a;
    }

    template <class Var0, class Var1, class Var2, class Var3>
    std::ostream& operator << (std::ostream& os, element_descriptor_impl<Var0,Var1,Var2,Var3> const& e) {
        os<<e.var0;
        os<<e.var1;
        os<<e.var2;
        os<<e.var3;
        return os;
    }


    namespace max_order_combined_helpers {

        template<unsigned int NK, unsigned int K>
        struct size_helper {
            static unsigned int const value = NK*size_helper<NK-1,K>::value;
        };

        template <unsigned int K>
        struct size_helper<K,K> {
            static unsigned int const value = 1;
        };

        template <unsigned int N>
        struct factorial {
            static unsigned int const value = N*factorial<N-1>::value;
        };

        template <>
        struct factorial<0> {
            static unsigned int const value = 1;
        };

        template <unsigned int N, unsigned int K>
        struct size {
            // N variables, max order K -> n+k-1 over k  = (n+k-1)! / ( (n-1)! k! ) combinations
            // Assuming N > 0
            static unsigned int const value = size_helper<N+K-1,K>::value/factorial<N-1>::value;
        };

    }


    template <class Variable, unsigned int Order>
    struct stride {
        static unsigned int const value = Order+1;
    };

    template <unsigned int Order>
    struct stride<no_variable,Order> {
        static unsigned int const value = 1;
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


    template <class Coeff, class OrderSpecification, class Var0, class Var1, class Var2, class Var3>
    struct storage {
    };
    
    
    // C - note for Andreas, sorry your element_descriptor costs too much, so I add an operator for direct access
    // C - I keep it for the print
    template <class Coeff, unsigned int Order, class Var0, class Var1, class Var2, class Var3>
    struct storage<Coeff, max_order_each<Order>, Var0, Var1, Var2, Var3> : public boost::array<Coeff, stride<Var0,Order>::value * stride<Var1,Order>::value * stride<Var2,Order>::value * stride<Var3,Order>::value> {
        typedef boost::array<Coeff, stride<Var0,Order>::value * stride<Var1,Order>::value * stride<Var2,Order>::value * stride<Var3,Order>::value>   base_type;
        typedef detail::element_descriptor_impl<Var0,Var1,Var2,Var3> element_descriptor;
//        typedef typename base_type::value_type  value_type;
        typedef typename base_type::size_type   size_type;
        // operator used for the inner product
        inline Coeff& operator ()(std::size_t i,std::size_t j,std::size_t k,std::size_t l) {
            return base_type::operator[](i*stride1*stride2*stride3 + j*stride2*stride3 + k*stride3 + l);            
        }

        inline Coeff const& operator ()(std::size_t i,std::size_t j,std::size_t k,std::size_t l) const{
            return base_type::operator[](i*stride1*stride2*stride3 + j*stride2*stride3 + k*stride3 + l);            
        }

        inline Coeff& operator()(element_descriptor const& e) {
            // TODO think of something smarter
            assert(e.var0.exp < stride0);
            assert(e.var1.exp < stride1);
            assert(e.var2.exp < stride2);
            assert(e.var3.exp < stride3);
            return base_type::operator[](e.var0.exp*stride1*stride2*stride3 + e.var1.exp*stride2*stride3 + e.var2.exp*stride3 + e.var3.exp);
        }
        
        inline Coeff const& operator()(element_descriptor const& e) const {
            assert(e.var0.exp < stride0);
            assert(e.var1.exp < stride1);
            assert(e.var2.exp < stride2);
            assert(e.var3.exp < stride3);
            return base_type::operator[](e.var0.exp*stride1*stride2*stride3 + e.var1.exp*stride2*stride3 + e.var2.exp*stride3 + e.var3.exp);
        }
 
      private:
        static unsigned int const stride0 = detail::stride<Var0,Order>::value; 
        static unsigned int const stride1 = detail::stride<Var1,Order>::value; 
        static unsigned int const stride2 = detail::stride<Var2,Order>::value; 
        static unsigned int const stride3 = detail::stride<Var3,Order>::value;
    };
    

    template <class Coeff, unsigned int Order, class Var0, class Var1, class Var2, class Var3>
    struct storage<Coeff, max_order_combined<Order>, Var0, Var1, Var2, Var3> : public boost::array<Coeff, max_order_combined_helpers::size<num_of_variables_helper<Var0,Var1,Var2,Var3>::value+1, Order>::value > {
        typedef boost::array<Coeff, max_order_combined_helpers::size<num_of_variables_helper<Var0,Var1,Var2,Var3>::value+1, Order>::value >    base_type;
        typedef detail::element_descriptor_impl<Var0, Var1, Var2, Var3> element_descriptor;
        typedef typename base_type::size_type   size_type;

        
        inline Coeff& operator ()(std::size_t i,std::size_t j,std::size_t k,std::size_t l) {
            return base_type::operator[](
                                         (i*(2*stride+3-i)*(2*stride*stride+6*stride+2 +i*i -2*stride*i - 3*i))/24 // Sum[1/6*(1*a^3+3*a^2+2*a)),{a,n-i+1,n}]
                                         + (j*( j*j - 3*j*(stride+1-i) + 3*(stride-i)*(stride+2-i) + 2))/6         // See <Var0,Var1,Var2>
                                         + (stride-i-j)*k - (k*k-k)/2 + l
            );            
        }
        
        inline Coeff const& operator ()(std::size_t i,std::size_t j,std::size_t k,std::size_t l) const{
            return base_type::operator[](
                                         (i*(2*stride+3-i)*(2*stride*stride+6*stride+2 +i*i -2*stride*i - 3*i))/24 // Sum[1/6*(1*a^3+3*a^2+2*a)),{a,n-i+1,n}]
                                         + (j*( j*j - 3*j*(stride+1-i) + 3*(stride-i)*(stride+2-i) + 2))/6         // See <Var0,Var1,Var2>
                                         + (stride-i-j)*k - (k*k-k)/2 + l
                                         );             
        }

        inline Coeff& operator()(element_descriptor const& e) {
            size_type const& i = e.var0.exp;
            size_type const& j = e.var1.exp;
            size_type const& k = e.var2.exp;
            size_type const& l = e.var3.exp;
            assert(i <= Order);
            assert(j <= Order-i);
            assert(k <= Order-i-j);
            assert(l <= Order-i-j-k);
            return base_type::operator[](
                  (i*(2*stride+3-i)*(2*stride*stride+6*stride+2 +i*i -2*stride*i - 3*i))/24 // Sum[1/6*(1*a^3+3*a^2+2*a)),{a,n-i+1,n}]
                + (j*( j*j - 3*j*(stride+1-i) + 3*(stride-i)*(stride+2-i) + 2))/6         // See <Var0,Var1,Var2>
                + (stride-i-j)*k - (k*k-k)/2 + l
                );
        }
        
        inline Coeff const& operator()(element_descriptor const& e) const {
            size_type const& i = e.var0.exp;
            size_type const& j = e.var1.exp;
            size_type const& k = e.var2.exp;
            size_type const& l = e.var3.exp;
            assert(i <= Order);
            assert(j <= Order-i);
            assert(k <= Order-i-j);
            assert(l <= Order-i-j-k);
            return base_type::operator[](
                  (i*(2*stride+3-i)*(2*stride*stride+6*stride+2 +i*i -2*stride*i - 3*i))/24
                + (j*( j*j - 3*j*(stride+1-i) + 3*(stride-i)*(stride+2-i) + 2))/6
                + (stride-i-j)*k - (k*k-k)/2 + l
                );
        }

      private:
        static unsigned int const stride = detail::stride<Var0,Order>::value; 
    };
    template <class Coeff, unsigned int Order, class Var0, class Var1, class Var2>
    struct storage<Coeff, max_order_combined<Order>, Var0, Var1, Var2, no_variable> : public boost::array<Coeff, max_order_combined_helpers::size<3+1,Order>::value> {
        typedef boost::array<Coeff, max_order_combined_helpers::size<3+1, Order>::value >    base_type;
        typedef detail::element_descriptor_impl<Var0, Var1, Var2, no_variable> element_descriptor;
        typedef typename base_type::size_type   size_type;

        
        inline Coeff& operator ()(std::size_t i,std::size_t j,std::size_t k,std::size_t l) {
            return base_type::operator[](
                                         (i*( i*i - 3*i*(stride+1) + 3*stride*(stride+2) + 2))/6  // Sum[(a^2+a)/2, {a,n-i+1,n}]
                                         + (stride-i)*j - (j*j-j)/2 + k                           // see <Var0,Var1>
                                         );            
        }
        
        inline Coeff const& operator ()(std::size_t i,std::size_t j,std::size_t k,std::size_t l) const{
            return base_type::operator[](
                                         (i*( i*i - 3*i*(stride+1) + 3*stride*(stride+2) + 2))/6  // Sum[(a^2+a)/2, {a,n-i+1,n}]
                                         + (stride-i)*j - (j*j-j)/2 + k                           // see <Var0,Var1>
                                         );             
        }
        
        inline Coeff& operator()(element_descriptor const& e) {
            size_type const& i = e.var0.exp;
            size_type const& j = e.var1.exp;
            size_type const& k = e.var2.exp;
            assert(i <= Order);
            assert(j <= Order-i);
            assert(k <= Order-i-j);
            return base_type::operator[](
                  (i*( i*i - 3*i*(stride+1) + 3*stride*(stride+2) + 2))/6  // Sum[(a^2+a)/2, {a,n-i+1,n}]
                + (stride-i)*j - (j*j-j)/2 + k                           // see <Var0,Var1>
                );
        }
        inline Coeff const& operator()(element_descriptor const& e) const {
            size_type const& i = e.var0.exp;
            size_type const& j = e.var1.exp;
            size_type const& k = e.var2.exp;
            assert(i <= stride);
            assert(j <= stride-i);
            assert(k <= stride-i-j);
            return base_type::operator[](
                  (i*( i*i - 3*i*(stride+1) + 3*stride*(stride+2) + 2))/6
                + (stride-i)*j - (j*j-j)/2 + k
                );
        }
         
      private:
        static unsigned int const stride = detail::stride<Var0,Order>::value; 
    };

    template <class Coeff, unsigned int Order, class Var0, class Var1>
    struct storage<Coeff, max_order_combined<Order>, Var0, Var1, no_variable, no_variable> : public boost::array<Coeff, max_order_combined_helpers::size<2+1,Order>::value> {
        typedef boost::array<Coeff, max_order_combined_helpers::size<2+1, Order>::value >    base_type;
        typedef detail::element_descriptor_impl<Var0, Var1, no_variable, no_variable> element_descriptor;
        typedef typename base_type::size_type   size_type;
        
        inline Coeff& operator ()(std::size_t i,std::size_t j,std::size_t k,std::size_t l) {
            return base_type::operator[](stride*i - (i*i - i)/2 + j);            
        }
        
        inline Coeff const& operator ()(std::size_t i,std::size_t j,std::size_t k,std::size_t l) const{
            return base_type::operator[](stride*i - (i*i - i)/2 + j);             
        }
        
        inline Coeff& operator()(element_descriptor const& e) {
            size_type const& i = e.var0.exp;
            size_type const& j = e.var1.exp;
            assert(i <= Order);
            assert(j <= Order-i);
            return base_type::operator[]( stride*i - (i*i - i)/2 + j); // Sum[a,{a,n-i+1,n}] + j
        }
        
        inline Coeff const& operator()(element_descriptor const& e) const {
            size_type const& i = e.var0.exp;
            size_type const& j = e.var1.exp;
            assert(i <= Order);
            assert(j <= Order-i);
            return base_type::operator[]( stride*i - (i*i - i)/2 + j);
        }
         
      private:
        static unsigned int const stride = detail::stride<Var0,Order>::value; 
    };
    
    template <class Coeff, unsigned int Order, class Var0>
    struct storage<Coeff, max_order_combined<Order>, Var0, no_variable, no_variable, no_variable> : public boost::array<Coeff, max_order_combined_helpers::size<1+1,Order>::value> {
        typedef boost::array<Coeff, max_order_combined_helpers::size<1+1, Order>::value >    base_type;
        typedef detail::element_descriptor_impl<Var0, no_variable, no_variable, no_variable> element_descriptor;
        typedef typename base_type::size_type   size_type;

        inline Coeff& operator ()(std::size_t i,std::size_t j,std::size_t k,std::size_t l) {
            return base_type::operator[](i);            
        }
        
        inline Coeff const& operator ()(std::size_t i,std::size_t j,std::size_t k,std::size_t l) const{
            return base_type::operator[](i);             
        }
        
        inline Coeff& operator()(element_descriptor const& e) {
            assert(e.var0.exp <= Order);
            return base_type::operator[](e.var0.exp);
        }
        inline Coeff const& operator()(element_descriptor const& e) const {
            assert(e.var0.exp <= Order);
            return base_type::operator[](e.var0.exp);
        }

    };
} // end namespace detail
} // end namespace vli
#endif //VLI_POLYNOMIAL_STORAGE_HPP
