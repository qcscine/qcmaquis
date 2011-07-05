/*
 *  monome_cpu.h
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 18.03.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */

#ifndef VLI_POLYNOME_CPU_H
#define VLI_POLYNOME_CPU_H
#include "boost/swap.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "detail/vli_polynomial_cpu_function_hooks.hpp"

#include "monome/monome.h"

#include <ostream>
#include <cmath>

namespace vli
{	    
    
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

// Andreas -> Should we renamed this class polynomial_cpu for convenience ? 
// As the class works with both implementation of VLI, keep presently polynomial
//
//C Sure. Rename it.
//C
//C I would write //C (for communication) if it is just a message to me
//C and not a real comment on the code.
//C That way we can clean the files easily afterwards.
//C
template<class Vli, int Order>
class polynomial{
public:
    typedef typename Vli::size_type size_type;      // Type of the exponents (has to be the same type as Vli::size_type)
    enum { max_order = Order};
    
    //   friend polynomial operator * <>(const polynomial& p, const monomial<Vli> & m);
    friend void poly_multiply <>(polynomial& result ,const polynomial& p1, const polynomial& p2);
    
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
        int n = memcmp((void*)&coeffs[0],(void*)&p.coeffs[0],Order*Order*Vli::size*sizeof(typename Vli::value_type));
        return (0 == n);
    }
    
    friend void swap(polynomial& p1, polynomial& p2)
    {
        boost::swap(p1.coeffs,p2.coeffs);
    }
    
    /**
     * Multiplies assign with coefficient
     * Mutliplies all elements the argument
     */
    polynomial& operator *= (Vli const& c)
    {
        for(int i=0; i<Order*Order;++i)
            coeffs[i] *= c;
        return *this;
    }
    
    
    polynomial& operator *= (monomial<Vli> const& c)
    {
        (*this) *= c.coeff;
        return *this;
    }
    
    /**
     * Access coefficient of monomial J^j_exp*h^h_exp
     */
    inline Vli operator ()(size_type j_exp, size_type h_exp) const
    {
        assert(j_exp < max_order);
        assert(h_exp < max_order);
        return coeffs[j_exp*max_order+h_exp];
    }
    
    /**
     * Access coefficient of monomial J^j_exp*h^h_exp
     */
    inline Vli& operator ()(size_type j_exp, size_type h_exp)
    {
        assert(j_exp < max_order);
        assert(h_exp < max_order);
        return coeffs[j_exp*max_order+h_exp];
    }
    
    void print(std::ostream& os) const
    {
        // TODO nice output
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

template<class Vli, int Order>//PUTAIN DE & !!!!!!! -> lol -> \o/
std::ostream& operator<<(std::ostream& os, polynomial<Vli, Order> const& p){
    p.print(os);
    return os;
}

} //end namespace
#endif
