#ifndef VLI_POLYNOMIAL_STORAGE_HPP
#define VLI_POLYNOMIAL_STORAGE_HPP

#include <vli/polynomial/polynomial_traits.hpp>
#include <vli/polynomial/detail/helpers.hpp>
#include <boost/array.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/contains.hpp>
#include <cassert>
#include <boost/static_assert.hpp>

namespace vli {
/* \cond I do not need this part in the doc*/
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
    inline typename boost::disable_if<boost::is_same<Var1,no_variable>,unsigned int&>::type exponent(element_descriptor_impl<Var0,Var1,Var2,Var3> & e, Var1 dummy) {
        return e.var1.exp;
    }
    
    template <class Var0, class Var1, class Var2, class Var3>
    inline typename boost::disable_if<boost::is_same<Var1,no_variable>,unsigned int const&>::type exponent(element_descriptor_impl<Var0,Var1,Var2,Var3> const& e, Var1 dummy) {
        return e.var1.exp;
    }
    
    template <class Var0, class Var1, class Var2, class Var3>
    inline typename boost::disable_if<boost::is_same<Var2,no_variable>,unsigned int&>::type exponent(element_descriptor_impl<Var0,Var1,Var2,Var3> & e, Var2 dummy) {
        return e.var2.exp;
    }
    
    template <class Var0, class Var1, class Var2, class Var3>
    inline typename boost::disable_if<boost::is_same<Var2,no_variable>,unsigned int const&>::type exponent(element_descriptor_impl<Var0,Var1,Var2,Var3> const& e, Var2 dummy) {
        return e.var2.exp;
    }
    
    template <class Var0, class Var1, class Var2, class Var3>
    inline typename boost::disable_if<boost::is_same<Var3,no_variable>,unsigned int&>::type exponent(element_descriptor_impl<Var0,Var1,Var2,Var3> & e, Var3 dummy) {
        return e.var3.exp;
    }

    template <class Var0, class Var1, class Var2, class Var3>
    inline typename boost::disable_if<boost::is_same<Var3,no_variable>,unsigned int const&>::type exponent(element_descriptor_impl<Var0,Var1,Var2,Var3> const& e, Var3 dummy) {
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


    template <class Coeff, class MaxOrder, int NumVars>
    struct storage;


    template <class Coeff, int Order, int NumVars>
    struct storage<Coeff, max_order_each<Order>, NumVars> : public boost::array<Coeff, num_coefficients<max_order_each<Order>,NumVars>::value> {
        typedef boost::array<Coeff, num_coefficients<max_order_each<Order>,NumVars>::value>   base_type;
        typedef typename base_type::size_type   size_type;
        inline Coeff& operator ()(std::size_t i, std::size_t j=0, std::size_t k=0, std::size_t l=0) {
            assert(i < stride0);
            assert(j < stride1);
            assert(k < stride2);
            assert(l < stride3);
            return base_type::operator[](i*stride1*stride2*stride3 + j*stride2*stride3 + k*stride3 + l);
        }

        inline Coeff const& operator ()(std::size_t i,std::size_t j=0, std::size_t k=0, std::size_t l=0) const {
            assert(i < stride0);
            assert(j < stride1);
            assert(k < stride2);
            assert(l < stride3);
            return base_type::operator[](i*stride1*stride2*stride3 + j*stride2*stride3 + k*stride3 + l);
        }

      private:
        static unsigned int const stride0 = detail::stride<0,NumVars,Order>::value;
        static unsigned int const stride1 = detail::stride<1,NumVars,Order>::value;
        static unsigned int const stride2 = detail::stride<2,NumVars,Order>::value;
        static unsigned int const stride3 = detail::stride<3,NumVars,Order>::value;
    };


    template <class Coeff, int Order>
    struct storage<Coeff, max_order_combined<Order>, 4> : public boost::array<Coeff, max_order_combined_helpers::size<4+1, Order>::value > {
        typedef boost::array<Coeff, max_order_combined_helpers::size<4+1, Order>::value >    base_type;
        typedef typename base_type::size_type   size_type;

        
        inline Coeff& operator ()(std::size_t i, std::size_t j=0, std::size_t k=0, std::size_t l=0) {
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
        
        inline Coeff const& operator ()(std::size_t i,std::size_t j=0, std::size_t k=0, std::size_t l=0) const {
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
      private:
        static unsigned int const stride = detail::stride<0,4,Order>::value;
    };
    template <class Coeff, int Order>
    struct storage<Coeff, max_order_combined<Order>, 3> : public boost::array<Coeff, max_order_combined_helpers::size<3+1, Order>::value> {
        typedef boost::array<Coeff, max_order_combined_helpers::size<3+1, Order>::value >    base_type;
        typedef typename base_type::size_type   size_type;

        
        inline Coeff& operator ()(std::size_t i, std::size_t j=0, std::size_t k=0, std::size_t l=0) {
            assert(i <= Order);
            assert(j <= Order-i);
            assert(k <= Order-i-j);
            return base_type::operator[](
                                         (i*( i*i - 3*i*(stride+1) + 3*stride*(stride+2) + 2))/6  // Sum[(a^2+a)/2, {a,n-i+1,n}]
                                         + (stride-i)*j - (j*j-j)/2 + k                           // see <Var0,Var1>
                                         );            
        }
        
        inline Coeff const& operator ()(std::size_t i, std::size_t j=0, std::size_t k=0, std::size_t l=0) const {
            assert(i <= Order);
            assert(j <= Order-i);
            assert(k <= Order-i-j);
            return base_type::operator[](
                                         (i*( i*i - 3*i*(stride+1) + 3*stride*(stride+2) + 2))/6  // Sum[(a^2+a)/2, {a,n-i+1,n}]
                                         + (stride-i)*j - (j*j-j)/2 + k                           // see <Var0,Var1>
                                         );             
        }
      private:
        static unsigned int const stride = detail::stride<0,3,Order>::value;
    };

    template <class Coeff, int Order>
    struct storage<Coeff, max_order_combined<Order>, 2> : public boost::array<Coeff, max_order_combined_helpers::size<2+1, Order>::value> {
        typedef boost::array<Coeff, max_order_combined_helpers::size<2+1, Order>::value >    base_type;
        typedef typename base_type::size_type   size_type;
        inline Coeff& operator ()(std::size_t i, std::size_t j=0, std::size_t k=0, std::size_t l=0) {
            assert(i <= Order);
            assert(j <= Order-i);
            return base_type::operator[](stride*i - (i*i - i)/2 + j);            
        }
        
        inline Coeff const& operator ()(std::size_t i, std::size_t j=0, std::size_t k=0, std::size_t l=0) const {
            assert(i <= Order);
            assert(j <= Order-i);
            return base_type::operator[](stride*i - (i*i - i)/2 + j);             
        }
      private:
        static unsigned int const stride = detail::stride<0,2,Order>::value;
    };
    
    template <class Coeff, int Order>
    struct storage<Coeff, max_order_combined<Order>, 1> : public boost::array<Coeff, max_order_combined_helpers::size<1+1, Order>::value> {
        typedef boost::array<Coeff, max_order_combined_helpers::size<1+1, Order>::value >    base_type;
        typedef typename base_type::size_type   size_type;

        inline Coeff& operator ()(std::size_t i, std::size_t j=0, std::size_t k=0, std::size_t l=0) {
            assert(i <= Order);
            return base_type::operator[](i);
        }
        
        inline Coeff const& operator ()(std::size_t i, std::size_t j=0, std::size_t k=0, std::size_t l=0) const {
            assert(i <= Order);
            return base_type::operator[](i);
        }
    };
} // end namespace detail
/* \endcond I do not need this part in the doc*/
} // end namespace vli
#endif //VLI_POLYNOMIAL_STORAGE_HPP
