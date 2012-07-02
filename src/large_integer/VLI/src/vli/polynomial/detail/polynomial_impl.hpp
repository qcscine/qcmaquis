
/***********************************************************************************
 * Copyright (C) 2012 by Andreas Hehn <hehn@phys.ethz.ch>                          *
 *                       ETH Zurich                                                *
 ***********************************************************************************/
#ifndef VLI_POLYNOMIAL_IMPL_HPP
#define VLI_POLYNOMIAL_IMPL_HPP

#include <vli/polynomial/detail/storage.hpp>
#include <vli/polynomial/polynomial_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <ostream>
#include <cassert>

#define POLYNOMIAL_CLASS polynomial<Coeff,OrderSpecification,Var0,Var1,Var2,Var3>

namespace vli {
namespace detail {
    namespace operations {
        struct plus_assign {
            template <class T, class T2>
            void operator()(T& t, T2 const& t2) {
                t += t2;
            }
        };

        struct minus_assign {
            template <class T, class T2>
            void operator()(T& t, T2 const& t2) {
                t -= t2;
            }
        };

        struct multiply_assign {
            template <class T, class T2>
            void operator()(T& t, T2 const& t2) {
                t *= t2;
            }
        };

        struct devide_assign {
            template <class T, class T2>
            void operator()(T& t, T2 const& t2) {
                t /= t2;
            }
        };
    } // end namespace operations 
    
    //
    // loop helper
    //
    template <class Polynomial>
    struct loop_helper {
    };
    
    template <class Coeff, unsigned int Order, class Var0, class Var1, class Var2, class Var3>
    struct loop_helper< polynomial<Coeff,max_order_each<Order>,Var0,Var1,Var2,Var3> > {
        typedef polynomial<Coeff,max_order_each<Order>,Var0,Var1,Var2,Var3> polynomial_type;
        typedef typename exponent_type<polynomial_type>::type               exponent_type;
        typedef typename element_descriptor<polynomial_type>::type          element_descriptor;

        template <class Operation>
        static void apply(polynomial_type const& p, Operation op) {
            for(exponent_type i=0; i < stride<Var0,Order>::value; ++i)
                for(exponent_type j=0; j < stride<Var1,Order>::value; ++j)
                    for(exponent_type k=0; k < stride<Var2,Order>::value; ++k)
                        for(exponent_type l=0; l < stride<Var3,Order>::value; ++l) 
                            op(p,element_descriptor(i,j,k,l));
        }
    };

    template <class Coeff, unsigned int Order, class Var0, class Var1, class Var2, class Var3>
    struct loop_helper< polynomial<Coeff,max_order_combined<Order>,Var0,Var1,Var2,Var3> > {
        typedef polynomial<Coeff,max_order_combined<Order>,Var0,Var1,Var2,Var3> polynomial_type;
        typedef typename exponent_type<polynomial_type>::type                   exponent_type;
        typedef typename element_descriptor<polynomial_type>::type              element_descriptor;

        template <class Operation>
        static void apply(polynomial_type const& p, Operation op) {
            element_descriptor el(0);
            apply_helper<Operation,0>(p, op, el, exponent_type(0), Var0() );
        }
        
        template <class Operation, unsigned int N, class VarN>
        static inline void apply_helper(polynomial_type const& p, Operation op, element_descriptor el, exponent_type order_sum, VarN d) {
            for(exponent_type i=0; i+order_sum < stride<VarN,Order>::value; ++i) {
                exponent(el,VarN()) = i;
                apply_helper<Operation,N+1>(p, op, el, order_sum+i, typename variable<polynomial_type,N+1>::type() );
            }
        }
        
        template <class Operation, unsigned int N>
        static inline void apply_helper(polynomial_type const& p, Operation op, element_descriptor const& el, exponent_type order_sum, no_variable d) {
            op(p,el);
        }
        
        template <class BinaryOperation, class ElementDescriptor>
        static void apply(polynomial_type& p, polynomial_type const& p2, BinaryOperation op, ElementDescriptor const& shift) {
            element_descriptor el(0);
            apply_helper<BinaryOperation,ElementDescriptor,0>(p, p2, op, shift, el, exponent_type(0), Var0() );
        }
        
        // These templates iterate over all variables starting from VarN
        template <class BinaryOperation, class ElementDescriptor, unsigned int N, class VarN>
        static inline void apply_helper(polynomial_type& p, polynomial_type const& p2, BinaryOperation op, ElementDescriptor const& shift, element_descriptor el, exponent_type order_sum, VarN d) {
            for(exponent_type i=exponent(shift,VarN()); i+order_sum < stride<VarN,Order>::value; ++i) {
                exponent(el,VarN()) = i;
                apply_helper<BinaryOperation,ElementDescriptor,N+1>(p, p2, op, shift, el, order_sum+i, typename variable<polynomial_type,N+1>::type() );
            }
        }
        
        template <class BinaryOperation, class ElementDescriptor, unsigned int N>
        static inline void apply_helper(polynomial_type& p, polynomial_type const& p2, BinaryOperation op, ElementDescriptor const& shift, element_descriptor const& el, exponent_type order_sum, no_variable d) {
            op(p,p2,shift,el);
        }
    };

    //
    // simple helper functions
    //
    template <class Coeff, class OrderSpecification, class Var0, class Var1, class Var2, class Var3>
    void init(POLYNOMIAL_CLASS & p, boost::false_type dummy) {
        // This is a non fundamental type -> it will be default constructed
    }
    
    template <class Coeff, class OrderSpecification, class Var0, class Var1, class Var2, class Var3>
    void init(POLYNOMIAL_CLASS & p, boost::true_type dummy) {
        // This is a fundamental type (e.g. unsigned int, double, ...) -> we have to initalize
        for(typename iterator<POLYNOMIAL_CLASS>::type it = p.begin(); it != p.end(); ++it)
            *it = Coeff();
    }
    
    template <class Coeff, class OrderSpecification, class Var0, class Var1, class Var2, class Var3, class Operation>
    void additive_op_assign(POLYNOMIAL_CLASS& p, POLYNOMIAL_CLASS const& p2, Operation op) {
        typename iterator<POLYNOMIAL_CLASS>::type it = p.begin();
        typename const_iterator<POLYNOMIAL_CLASS>::type it2 = p2.begin();
        while(it != p.end())
            op(*it++, *it2++);
    }

    template <class Coeff, class OrderSpecification, class Var0, class Var1, class Var2, class Var3, class MCoeff, class MVar0, class MVar1, class MVar2, class MVar3, class Operation>
    void additive_op_assign(POLYNOMIAL_CLASS& p, monomial<MCoeff,MVar0,MVar1,MVar2,MVar3> const& m, Operation op) {
        op(p(typename element_descriptor<POLYNOMIAL_CLASS>::type(m)), m.c_);
    }
    
    template <class Coeff, class OrderSpecification, class Var0, class Var1, class Var2, class Var3, class Operation>
    void additive_op_assign(POLYNOMIAL_CLASS& p, Coeff const& c, Operation op) {
        op(*p.begin(),c);
    }

    template <class Coeff, class OrderSpecification, class Var0, class Var1, class Var2, class Var3, class Operation>
    void additive_op_assign(POLYNOMIAL_CLASS& p, typename boost::enable_if<boost::is_same<Coeff,int>,int>::type a, Operation op) {
        op(*p.begin(),a);
    }
    
    template <class Coeff, class OrderSpecification, class Var0, class Var1, class Var2, class Var3, class Operation, class Coeff2>
    void multiplicative_op_assign(POLYNOMIAL_CLASS& p, Coeff2 const& c, Operation op) {
        for(typename iterator<POLYNOMIAL_CLASS>::type it=p.begin(); it != p.end(); ++it)
            op(*it, c);
    }
    
//    template <class Coeff, class OrderSpecification, class Var0, class Var1, class Var2, class Var3, class Operation>
//    void multiplicative_op_assign(POLYNOMIAL_CLASS& p, typename boost::enable_if<boost::is_same<Coeff,int>,int>::type a, Operation op) {
//        for(typename iterator<POLYNOMIAL_CLASS>::type it=p.begin(); it != p.end(); ++it)
//            op(*it, a);
//    }
    
    template <class Polynomial>
    struct print_helper {
        struct print_element {
            print_element(std::ostream& os)
            : os_(os) {
            }
            
            void operator()(Polynomial const& p, typename element_descriptor<Polynomial>::type const& e) {
                if( !(p(e) < typename Polynomial::value_type(0)) ) 
                    os_<<"+";
                os_<<p(e)<<e<<" " << std::endl;
            }
            
          private:
            std::ostream&       os_;
        };
        static void apply(std::ostream& os, Polynomial const& p) {
            loop_helper<Polynomial>::apply(p,print_element(os));
        }
    };
    
    template <class T>
    struct multiply_factor_op {
        multiply_factor_op(T const& factor)
        : c_(factor){
        }
        template <class Polynomial, class ElementDescriptor>
        inline void operator()(Polynomial& p, Polynomial const& p2, ElementDescriptor const& shift, typename element_descriptor<Polynomial>::type const& element) {
            p(element) = p2(element-shift) * c_;
        }
        private:
        T c_;
    };
  
    template <class Coeff, class OrderSpecification, class Var0, class Var1, class Var2, class Var3>
    void negate(POLYNOMIAL_CLASS& p) {
        for(typename iterator<POLYNOMIAL_CLASS>::type it=p.begin(); it != p.end(); ++it)
            negate_inplace(*it);
    }
    
    template <class Coeff, class OrderSpecification, class Var0, class Var1, class Var2, class Var3>
    bool is_zero_helper(POLYNOMIAL_CLASS const& p) {
        for(typename const_iterator<POLYNOMIAL_CLASS>::type it=p.begin(); it != p.end(); ++it)
            if (!is_zero(*it))
                return false;
        return true;
    }
    
    template <class Polynomial>
    struct equal_helper {
        bool operator()(Polynomial const& p1, Polynomial const& p2) {
            typedef typename const_iterator<Polynomial>::type const_iterator;
            const_iterator it(p1.begin()), it2(p2.begin());
            while(it != p1.end()) {
                if (!(*it == *it2))
                    return false;
                ++it;
                ++it2;
            }
            return true;
        }
    };

    

    //
    // monomial multiply helper
    //
    template <class Polynomial>
    struct multiply_assign_monomial_helper{
    };
   
    template <class Coeff, unsigned int Order, class Var0, class Var1, class Var2, class Var3>
    struct multiply_assign_monomial_helper< polynomial<Coeff,max_order_each<Order>,Var0,Var1,Var2,Var3> > {
        typedef polynomial<Coeff,max_order_each<Order>,Var0,Var1,Var2,Var3> polynomial_type;
        typedef typename exponent_type<polynomial_type>::type               exponent_type;

        template <class MCoeff, class MVar0, class MVar1, class MVar2, class MVar3>
        static void apply(polynomial_type& p, monomial<MCoeff,MVar0,MVar1,MVar2,MVar3> const& m) {
           //TODO optimize: this operation should be possible without creating a new copy (is it worth the effort?)
           polynomial_type p2;
           for(exponent_type i=exponent(m,Var0()); i < stride<Var0,Order>::value; ++i)
               for(exponent_type j=exponent(m,Var1()); j < stride<Var1,Order>::value; ++j)
                   for(exponent_type k=exponent(m,Var2()); k < stride<Var2,Order>::value; ++k)
                       for(exponent_type l=exponent(m,Var3()); l < stride<Var2,Order>::value; ++l) 
                           p2(i,j,k,l) = m.c_ * p(i-exponent(m,Var0()),j-exponent(m,Var1()),k-exponent(m,Var2()),l-exponent(m,Var3()));
           swap(p2,p);
        }
    };
    
    template <class Coeff, unsigned int Order, class Var0>
    struct multiply_assign_monomial_helper< polynomial<Coeff, max_order_each<Order>, Var0, no_variable, no_variable, no_variable> > {
        typedef polynomial<Coeff, max_order_each<Order>, Var0, no_variable, no_variable, no_variable>   polynomial_type;
        typedef typename polynomial_type::reverse_iterator                                              reverse_iterator;
        typedef typename polynomial_type::value_type                                                    value_type;
        
        template <class MCoeff>
        static void apply(polynomial_type& p, monomial<MCoeff, Var0, no_variable, no_variable, no_variable> const& m) {
           for(reverse_iterator it = p.rbegin() + exponent(m,Var0()); it !=p.rend(); ++it)
               *(it-m.var0.exp) = m.c_ * (*it);
           for(reverse_iterator it = p.rend() - exponent(m,Var0()); it != p.rend(); ++it)
               *(it) = value_type(0);
        }
    };
    
    template <class Coeff, unsigned int Order, class Var0, class Var1, class Var2, class Var3>
    struct multiply_assign_monomial_helper< polynomial<Coeff,max_order_combined<Order>,Var0,Var1,Var2,Var3> > {
        typedef polynomial<Coeff,max_order_combined<Order>,Var0,Var1,Var2,Var3> polynomial_type;
        typedef typename exponent_type<polynomial_type>::type               exponent_type;
        template <class MVar0, class MVar1, class MVar2, class MVar3>
        static void apply(polynomial_type& p, monomial<Coeff,MVar0,MVar1,MVar2,MVar3> const& m) {
           polynomial_type p2;
           loop_helper<polynomial_type>::apply(p2,p,multiply_factor_op<Coeff>(m.c_),m);
           swap(p2,p); 
        }
    };
    
    //
    // polynomial_multiply_helper
    //
    template <class Polynomial>
    struct polynomial_multiply_helper {
    };
    
 

    template <class Coeff, unsigned int Order, class Var0, class Var1, class Var2, class Var3>
    struct polynomial_multiply_helper< polynomial<Coeff,max_order_each<Order>,Var0,Var1,Var2,Var3> > {
        typedef polynomial<Coeff,max_order_each<Order>,Var0,Var1,Var2,Var3>     polynomial_type;
        typedef typename polynomial_multiply_result_type<polynomial_type>::type result_type;
        typedef typename element_descriptor<result_type>::type                  result_element_descriptor;
        typedef typename exponent_type<polynomial_type>::type                   exponent_type;
        typedef typename element_descriptor<polynomial_type>::type              element_descriptor;
        static result_type apply(polynomial_type const& p1 , polynomial_type const& p2) {
            // TODO optimize
            result_type result;

            for(exponent_type i = 0; i < stride<Var0,Order>::value; ++i)
                for(exponent_type i2 = 0; i2 < stride<Var0,Order>::value; ++i2)
                    for(exponent_type j = 0; j < stride<Var1,Order>::value; ++j)
                        for(exponent_type j2 = 0; j2 < stride<Var1,Order>::value; ++j2)                        
                            for(exponent_type k = 0; k < stride<Var2,Order>::value; ++k)
                                for(exponent_type k2 = 0; k2 < stride<Var2,Order>::value; ++k2)
                                    for(exponent_type l = 0; l < stride<Var3,Order>::value; ++l)
                                        for(exponent_type l2 = 0; l2 < stride<Var3,Order>::value; ++l2)
                                            multiply_add(result(i+i2,j+j2,k+k2,l+l2), p1(i,j,k,l), p2(i2,j2,k2,l2));
            

             return result;
        };
    };

    template <class Coeff, unsigned int Order, class Var0, class Var1, class Var2, class Var3>
    struct polynomial_multiply_helper< polynomial<Coeff,max_order_combined<Order>,Var0,Var1,Var2,Var3> > {
        typedef polynomial<Coeff,max_order_combined<Order>,Var0,Var1,Var2,Var3> polynomial_type;
        typedef typename polynomial_multiply_result_type<polynomial_type>::type result_type;
        typedef typename element_descriptor<result_type>::type                  result_element_descriptor;
        typedef typename exponent_type<polynomial_type>::type                   exponent_type;
        typedef typename element_descriptor<polynomial_type>::type              element_descriptor;
        static result_type apply(polynomial_type const& p1 , polynomial_type const& p2) {
            result_type result;
            apply_helper<0>(result,p1,p2,result_element_descriptor(0),element_descriptor(0),element_descriptor(0),0,0,Var0());
            return result;
        };

        // TODO this is a naive loop algorithm using element-wise access -> optimize
        template <unsigned int N,class VarN>
        static void apply_helper(result_type& result, polynomial_type const& p1, polynomial_type const& p2, result_element_descriptor er, element_descriptor e1, element_descriptor e2, exponent_type order_sum, exponent_type order_sum2, VarN d) {
            for(exponent_type i=0; i+order_sum < stride<VarN,Order>::value; ++i) {
                exponent(e1,VarN()) = i;
                for(exponent_type i2=0; i2+order_sum2 < stride<VarN,Order>::value; ++i2) {
                    exponent(e2,VarN()) = i2;
                    exponent(er,VarN()) = i+i2;
                    apply_helper<N+1>(result,p1,p2,er,e1,e2,order_sum+i,order_sum2+i2,typename variable<polynomial_type,N+1>::type());
                }
            }
        }

        template <unsigned int N>
        static inline void apply_helper(result_type& result, polynomial_type const& p1, polynomial_type const& p2, result_element_descriptor const& er, element_descriptor const& e1, element_descriptor const& e2, exponent_type order_sum, exponent_type order_sum2, no_variable d) {
            multiply_add(result(er), p1(e1), p2(e2));
        }
    };
    
    //
    // polynomial_multiply_helper
    //
    template <class Polynomial>
    struct polynomial_multiply_keep_order_helper {
    };
    
    template <class Coeff, unsigned int Order, class Var0, class Var1, class Var2, class Var3>
    struct polynomial_multiply_keep_order_helper< polynomial<Coeff,max_order_each<Order>,Var0,Var1,Var2,Var3> > {
        typedef polynomial<Coeff,max_order_each<Order>,Var0,Var1,Var2,Var3>                 polynomial_type;
        typedef typename polynomial_multiply_keep_order_result_type<polynomial_type>::type  result_type;
        typedef typename element_descriptor<result_type>::type                              result_element_descriptor;
        typedef typename exponent_type<polynomial_type>::type                               exponent_type;
        typedef typename element_descriptor<polynomial_type>::type                          element_descriptor;
        static typename polynomial_multiply_keep_order_result_type<polynomial_type>::type apply(polynomial_type const& p1, polynomial_type const& p2) {
            result_type result;
            for(exponent_type i = 0; i < stride<Var0,Order>::value; ++i)
                for(exponent_type i2 = 0; i2+i < stride<Var0,Order>::value; ++i2)
                    for(exponent_type j = 0; j < stride<Var1,Order>::value; ++j)
                        for(exponent_type j2 = 0; j2+j < stride<Var1,Order>::value; ++j2)
                            for(exponent_type k = 0; k < stride<Var2,Order>::value; ++k)
                                for(exponent_type k2 = 0; k2+k < stride<Var2,Order>::value; ++k2)
                                    for(exponent_type l = 0; l < stride<Var3,Order>::value; ++l)
                                        for(exponent_type l2 = 0; l2+l < stride<Var3,Order>::value; ++l2)
                                              multiply_add(result(i+i2,j+j2,k+k2,l+l2), p1(i,j,k,l), p2(i2,j2,k2,l2));
            return result;
        }
    };
    
    template <class Coeff, unsigned int Order, class Var0, class Var1, class Var2, class Var3>
    struct polynomial_multiply_keep_order_helper< polynomial<Coeff,max_order_combined<Order>,Var0,Var1,Var2,Var3> > {
        typedef polynomial<Coeff,max_order_combined<Order>,Var0,Var1,Var2,Var3> polynomial_type;
        typedef typename polynomial_multiply_keep_order_result_type<polynomial_type>::type result_type;
        typedef typename element_descriptor<result_type>::type                  result_element_descriptor;
        typedef typename exponent_type<polynomial_type>::type                   exponent_type;
        typedef typename element_descriptor<polynomial_type>::type              element_descriptor;
        static result_type apply(polynomial_type const& p1 , polynomial_type const& p2) {
            result_type result;
            apply_helper<0>(result,p1,p2,result_element_descriptor(0),element_descriptor(0),element_descriptor(0),0,0,Var0());
            return result;
        };
// TODO correct implementation
        // TODO this is a naive loop algorithm using element-wise access -> optimize
        template <unsigned int N,class VarN>
        static void apply_helper(result_type& result, polynomial_type const& p1, polynomial_type const& p2, result_element_descriptor er, element_descriptor e1, element_descriptor e2, exponent_type order_sum, exponent_type order_sum2, VarN d) {
            for(exponent_type i=0; i+order_sum < stride<VarN,Order>::value; ++i) {
                exponent(e1,VarN()) = i;
                for(exponent_type i2=0; i2+order_sum2 < stride<VarN,Order>::value; ++i2) {
                    exponent(e2,VarN()) = i2;
                    exponent(er,VarN()) = i+i2;
                    apply_helper<N+1>(result,p1,p2,er,e1,e2,order_sum+i,order_sum2+i2,typename variable<polynomial_type,N+1>::type());
                }
            }
        }

        template <unsigned int N>
        static inline void apply_helper(result_type& result, polynomial_type const& p1, polynomial_type const& p2, result_element_descriptor const& er, element_descriptor const& e1, element_descriptor const& e2, exponent_type order_sum, exponent_type order_sum2, no_variable d) {
            multiply_add(result(er), p1(e1), p2(e2));
        }
    };

 
} // end namespace detail
} // end namespace vli

#undef POLYNOMIAL_CLASS

#endif //VLI_POLYNOMIAL_IMPL_HPP
