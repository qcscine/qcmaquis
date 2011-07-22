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
#include "function_hooks/vli_polynomial_cpu_function_hooks.hpp"

#include "monome/monome.h"

#include <ostream>
#include <cmath>

namespace vli
{	    
    
template<class Vli, int Order>
class polynomial_cpu;

/**
 * Multiplication of two polynomial_cpus
 */
template<class Vli, int Order>
polynomial_cpu<Vli, Order> operator * (polynomial_cpu<Vli, Order> const& p1, polynomial_cpu<Vli, Order> const& p2)
{
    polynomial_cpu<Vli, Order> result;
    poly_multiply(result, p1, p2);
    return result;
}


template<class Vli, int Order>
class polynomial_cpu{
public:
    typedef typename Vli::size_type size_type;      // Type of the exponents (has to be the same type as Vli::size_type)
    typedef typename Vli::value_type value_type;
    enum { size = Vli::size};
    enum { max_order = Order};
    
    //   friend polynomial_cpu operator * <>(const polynomial_cpu& p, const monomial<Vli> & m);
    friend void poly_multiply <>(polynomial_cpu& result ,const polynomial_cpu& p1, const polynomial_cpu& p2);
    
    polynomial_cpu(){
        for(int i=0; i<Order*Order;++i)
            coeffs[i]=Vli();
    }
    
    polynomial_cpu(const polynomial_cpu& p){
        for(int i=0; i<Order*Order;++i)
            coeffs[i]=p.coeffs[i];
    }
    
    polynomial_cpu& operator = (polynomial_cpu p)
    {
        swap(*this,p);
        return *this;
    }
    
    /**
     * Plus assign with a polynomial_cpu
     */
    polynomial_cpu& operator += (polynomial_cpu const& p)
    {
        for(int i=0; i<Order*Order;++i)
            coeffs[i]+=p.coeffs[i];
        return *this;
    }
    
    template <typename T>
    polynomial_cpu& operator += (T const& t)
    {
        coeffs[0]+=t;
        return *this;
    }
    
    /**
     * Minus assign with a polynomial_cpu
     */
    polynomial_cpu& operator -= (polynomial_cpu const& p)
    {
        for(int i=0; i<Order*Order;++i)
            coeffs[i]-=p.coeffs[i];
        
        return *this;
    }

    bool operator==(polynomial_cpu const& p) const
    {
        int n = memcmp((void*)&coeffs[0],(void*)&p.coeffs[0],Order*Order*Vli::size*sizeof(typename Vli::value_type));
        return (0 == n);
    }
    
    friend void swap(polynomial_cpu& p1, polynomial_cpu& p2)
    {
        boost::swap(p1.coeffs,p2.coeffs);
    }
    
    /**
     * Multiplies assign with coefficient
     * Mutliplies all elements the argument
     */
    polynomial_cpu& operator *= (Vli const& c)
    {
        for(int i=0; i<Order*Order;++i)
            coeffs[i] *= c;
        return *this;
    }
    
    
    polynomial_cpu& operator *= (monomial<Vli> const& c)
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
        for(std::size_t i = 0; i < Order ; i++){
            for(std::size_t j = 0; j < Order ; j++){
                os << "Coeff (j,h) = " << i <<" "<<j<< std::endl;
                os <<coeffs[i*Order+j] << std::endl;
            }
        }
    }
     
private:
    Vli coeffs[Order*Order];
};

/*
 * Multiplication of a monomial with a polynomial_cpu
 */
template<class Vli, int Order>
polynomial_cpu<Vli, Order> operator * (monomial<Vli> const& m,polynomial_cpu<Vli, Order> const& p)
{
    return p * m;
}

/**
 * Multiplication of a polynomial_cpu with a factor
 */
template<class Vli, int Order>
polynomial_cpu<Vli, Order> operator * (polynomial_cpu<Vli, Order> p, Vli const& c)
{
    p *= c;
    return p;
}

template<class Vli, int Order>
polynomial_cpu<Vli, Order> operator * (Vli const& c, polynomial_cpu<Vli, Order> const& p)
{
    return p * c;
}

template<class Vli, int Order> 
std::ostream& operator<<(std::ostream& os, polynomial_cpu<Vli, Order> const& p){
    p.print(os);
    return os;
}

} //end namespace
#endif
