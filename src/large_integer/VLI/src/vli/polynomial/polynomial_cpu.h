#ifndef VLI_POLYNOME_CPU_H
#define VLI_POLYNOME_CPU_H
#include "vli/function_hooks/vli_polynomial_cpu_function_hooks.hpp"
#include "vli/polynomial/monomial.hpp"
#include "vli/vli_cpu.h"

#include <boost/swap.hpp>
#include <ostream>
#include <cmath>
#include <cassert>

namespace vli
{	    

template<class Vli, unsigned int Order>          
void swap(polynomial_cpu<Vli,Order>& p1, polynomial_cpu<Vli,Order>& p2){
    using boost::swap;
    swap(p1.coeffs_,p2.coeffs_);
}

template<class Vli, unsigned int Order>
class polynomial_cpu{
public:
    typedef unsigned int exponent_type;      // Type of the exponents (has to be the same type as Vli::size_type)
    typedef Vli value_type;
    enum { max_order = Order};
        
    polynomial_cpu();
    explicit polynomial_cpu(int i);
    polynomial_cpu(const polynomial_cpu& p);    
    polynomial_cpu& operator = (polynomial_cpu p);
    
    polynomial_cpu& operator += (polynomial_cpu const& p);    
    polynomial_cpu& operator += (monomial<Vli> const& m);
    polynomial_cpu& operator += (int a);

    polynomial_cpu& operator -= (polynomial_cpu const& p);    
    polynomial_cpu& operator -= (monomial<Vli> const& m);
    polynomial_cpu& operator -= (int a);
    
    polynomial_cpu operator - () const;
    void negate();

    bool operator==(polynomial_cpu const& p) const;
    bool is_zero() const;
    
    friend void swap<>(polynomial_cpu<Vli,Order>& p1, polynomial_cpu<Vli,Order>& p2);

    polynomial_cpu& operator *= (monomial<Vli> const& m);
    polynomial_cpu& operator *= (Vli const& c);
    polynomial_cpu& operator *= (int c);
    
    inline Vli const& operator ()(exponent_type j_exp, exponent_type h_exp) const;
    inline Vli& operator ()(exponent_type j_exp, exponent_type h_exp);
    
    void print(std::ostream& os) const;
     
    Vli coeffs_[Order*Order];
};

template <class Vli, unsigned int Order> 
void poly_multiply (polynomial_cpu<Vli,2*Order>& result, const polynomial_cpu<Vli,Order>& p1, const polynomial_cpu<Vli,Order>& p2);   


template<class T, class Vli, unsigned int Order>
polynomial_cpu<Vli, Order> operator * (monomial<T> const& m,polynomial_cpu<Vli, Order> const& p);

template<class Vli, unsigned int Order>
polynomial_cpu<Vli, Order> operator * (polynomial_cpu<Vli, Order> p, Vli const& c);

template<class Vli, unsigned int Order>
polynomial_cpu<Vli, Order> operator * (polynomial_cpu<Vli, Order> p, int c);

template<class Vli, unsigned int Order>
polynomial_cpu<Vli, Order> operator * (Vli const& c, polynomial_cpu<Vli, Order> const& p);

template<class Vli, unsigned int Order>
polynomial_cpu<Vli, Order> operator * (int c, polynomial_cpu<Vli, Order> const& p);

template <class BaseInt, std::size_t Size, unsigned int Order> // C - the VliOut is twice larger
polynomial_cpu<vli_cpu<BaseInt, 2*Size>, 2*Order> operator * (polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const& p1, polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const& p2);
    
template <class BaseInt, std::size_t Size, unsigned int Order, class T>
polynomial_cpu<vli_cpu<BaseInt, Size>, Order> operator * (polynomial_cpu<vli_cpu<BaseInt, Size>, Order> p, monomial<T> const& m);

template<class Vli, unsigned int Order> 
std::ostream& operator<<(std::ostream& os, polynomial_cpu<Vli, Order> const& p);
} //end namespace
#include "vli/polynomial/polynomial_cpu.hpp"

#endif
