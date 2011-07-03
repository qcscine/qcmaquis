/*
 *  monome.h
 *  vli
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 18.03.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */

#ifndef VLI_MONOME_H
#define VLI_MONOME_H
#include "detail/vli_polynomial_cpu_function_hooks.hpp"
//#include "detail/vli_polynomial_gpu_function_hooks.hpp"
#include "boost/swap.hpp"

#include <ostream>
#include <cmath>


#define POLYNOMIAL_MAX_ORDER 4 // for testing

namespace vli
{	
	template <class Vli>
    struct monomial
    {
        enum { vli_size = Vli::vli_size };
        typedef typename Vli::size_type size_type;	
        typedef typename Vli::value_type value_type;
        /**
         * Constructor: Creates a monomial 1*J^j_exp*h^h_exp
         */
        explicit monomial(size_type j_exp = 0, size_type h_exp = 0)
        :j_exp(j_exp), h_exp(h_exp){
        }
        
        monomial& operator *= (Vli const& c);
        monomial& operator *= (value_type c);

        value_type& operator[](size_type i){
            return coeff[i];
        }
        
        size_type j_exp;
        size_type h_exp;
        Vli coeff;                
    };
    
    template<class Vli, int Order>
	class polynomial;
    
    /**
     * Multiplication of two polynomials
     */
   	template<class Vli, int Order>
    polynomial<Vli, Order> operator * (polynomial<Vli, Order> const& p1, polynomial<Vli, Order> const& p2)
    {
        polynomial<Vli, Order> result;
        poly_multiply(result, p1, p2);
        return result;
    }
    
	template<class Vli, int Order>
	class polynomial{
	public:
     //   friend polynomial operator * <>(const polynomial& p, const monomial<Vli> & m);
        
        friend polynomial operator * <>(const polynomial& p, const polynomial& m);
        friend void poly_multiply <>(polynomial& result ,const polynomial& p1, const polynomial& p2);
		
        typedef typename Vli::size_type size_type;	
        typedef typename Vli::value_type value_type;
        enum { vli_size  = Vli::vli_size };
        enum { max_order = Order};
        
        polynomial(){
            for(int i=0; i<Order*Order;++i)
                coeffs[i]=Vli();
        }

        polynomial(const polynomial& p){
            for(int i=0; i<Order*Order;++i)
                coeffs[i]=p.coeffs[i];
        }
        
        polynomial& operator = (polynomial p)
        {
            swap(*this,p);
            return *this;
        }
  
        /**
         * Plus assign with a polynomial
         */
        polynomial& operator += (polynomial const& p)
        {
            for(int i=0; i<Order*Order;++i)
                coeffs[i]+=p.coeffs[i];
            
            return *this;
        }
        
        bool operator==(polynomial const& p) const
        {
            int n = memcmp((void*)&coeffs[0],(void*)&p.coeffs[0],Order*Order*vli_size*sizeof(value_type));
			return (0 == n);
        }
        
        friend void swap(polynomial& p1, polynomial& p2)
        {
            boost::swap(p1.coeffs,p2.coeffs);
        }
        
        /**
         * Access coefficient of monomial J^j_exp*h^h_exp
         */
        inline Vli operator ()(unsigned int j_exp, unsigned int h_exp) const
        {
            assert(j_exp < max_order);
            assert(h_exp < max_order);
            return coeffs[j_exp*max_order+h_exp];
        }
        
        /**
         * Access coefficient of monomial J^j_exp*h^h_exp
         */
        inline Vli& operator ()(unsigned int j_exp, unsigned int h_exp)
        {
            assert(j_exp < max_order);
            assert(h_exp < max_order);
            return coeffs[j_exp*max_order+h_exp];
        }
        
        void print(std::ostream& os) const
        {
            for(int i=0; i<Order*Order;++i)
                std::cout<<coeffs[i] << " ";
		}
        
        private:
        Vli coeffs[Order*Order];
    };
    
    
    
    /*
     * Multiplication of a monomial with a polynomial
     */
   	template<class Vli, int Order>
    polynomial<Vli, Order> operator * (monomial<Vli> const& m,polynomial<Vli, Order> const& p)
    {
        return p * m;
    }
    
    /**
     * Multiplication of a polynomial with a factor
     */
    template<class Vli, int Order>
    polynomial<Vli, Order> operator * (polynomial<Vli, Order> p, Vli const& c)
    {
        p *= c;
        return p;
    }
    
    template<class Vli, int Order>
    polynomial<Vli, Order> operator * (Vli const& c, polynomial<Vli, Order> const& p)
    {
        return p * c;
    }
    
    template<class Vli, int Order>//PUTAIN DE & !!!!!!!
    std::ostream& operator<<(std::ostream& os, polynomial<Vli, Order> const& p){
        p.print(os);
        return os;
    }
    
    
}

#endif //VLI_MONOME_H
