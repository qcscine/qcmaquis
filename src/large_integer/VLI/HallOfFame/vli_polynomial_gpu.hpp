z/*****************************************************************************
 *
 * A minimal implementation of a polynomial in two symbolic variables 'J','h'.
 *
 * (C) 2010 by Andreas Hehn (hehn@phys.ethz.ch)
 *
 *****************************************************************************/

#ifndef HP2C__POLYNOMIAL
#define HP2C__POLYNOMIAL

#include <vector>
#include <iostream>
#include "vli_gpu/vli_polynomial_gpu.hpp"


namespace vli
{

/**
  * A monomial class
  */

template <class VLI>
struct monomial
{
    typedef typname VLI::size_type size_type j_exp;
    typedef typname VLI::size_type size_type h_exp;
    VLI coeff_;

    /**
      * Constructor: Creates a monomial 1*J^j_exp*h^h_exp
      */
    explicit monomial(size_type j_exp = 0, size_type h_exp = 0)
        : j_exp(j_exp), h_exp(h_exp), coeff_(1)
    {
    }

    monomial& operator *= (VLI const& c)
    {
        coeff *= c;
        return *this;
    }

    monomial& operator *= (int c)
    {
        coeff *= c;
        return *this;
    }
};

template <typename CoeffType, typename T>
monomial<CoeffType> operator * (monomial<CoeffType> m, T const& t)
{
    m *= t;
    return m;
}

template <typename CoeffType, typename T>
monomial<CoeffType> operator * (T const& t, monomial<CoeffType> const& m)
{
    return m*t;
}



template <class VLI>
struct polynomial
{
    typedef typename VLI::value_type value_type
    
    vli_vector_gpu<value_type> coeffs_; 
    
    
    enum {max_order = POLYNOMIAL_MAX_ORDER};

    friend polynomial operator *<> (polynomial const&, polynomial const&);

        /**
          * Constructor: Creates a polynomial
          * p(J,h) = c_0 + c_{1,0} * J^1 * h^0 + c_{2,0} * J^2 * h^0 + ...
          *             + c_{max_order-1,max_order-1} * (J*h)^(max_order-1)
          * where all c_* are of type CoeffType and set to 0.
          */
    polynomial():coeffs_(max_order*max_order){}

    polynomial(polynomial const& p):coeffs_(p.coeffs_){}

    polynomial& operator = (polynomial p){
        swap(p);
        return *this;
    }

    void swap(polynomial& p){
        coeffs.swap(p.coeffs);
    }

    void print(std::ostream& o) const{
        for(std::size_t je = 0; je < max_order; ++je)
            for(std::size_t he = 0; he < max_order; ++he)
                if( coeffs[je*max_order+he] != CoeffType(0))
                    o<<"+"<<coeffs[je*max_order+he]<<"*J^"<<je<<"*h^"<<he;
    }

    polynomial& operator += (polynomial const& p){
        this->coeffs_+=p.coeffs_;
        return *this;
    }

    polynomial& operator += (VLI const& c){
        coeffs[0] += c;
        return *this;
    }

    polynomial& operator += (monomial<VLI> const& m){
        assert(m.j_exp < max_order);
        assert(m.h_exp < max_order);
        coeffs[m.j_exp*max_order+m.h_exp] += m.coeff;
        return *this;
    }

    polynomial& operator *= (VLI const& c){
        this->coeffs_+=p.coeffs_;
        return *this;
    }

    /**
    * Access coefficient of monomial J^j_exp*h^h_exp
    */

    inline CoeffType operator ()(unsigned int j_exp, unsigned int h_exp) const{
        assert(j_exp < max_order);
        assert(h_exp < max_order);
        return coeffs[j_exp*max_order+h_exp];
    }
        
    /**
    * Access coefficient of monomial J^j_exp*h^h_exp
    */
    inline CoeffType& operator ()(unsigned int j_exp, unsigned int h_exp){
        assert(j_exp < max_order);
        assert(h_exp < max_order);
        return coeffs[j_exp*max_order+h_exp];
    }
};

template <class VLI>
std::ostream& operator <<(std::ostream& o, polynomial<VLI> const& p)
{
    p.print(o);
    return o;
}

template <class VLI>
polynomial<VLI> operator * (monomial<VLI> const& m,polynomial<VLI> const& p)
{
    return p * m;
}

template <class VLI>
polynomial<VLI> operator * (polynomial<VLI> p, CoeffType const& c)
{
    p *= c;
    return p;
}

template <class VLI>
polynomial<VLI> operator * (CoeffType const& c, polynomial<VLI> const& p)
{
    return p * c;
}

}



//
// *****************************************************************************
//

/**
 * A polynomial class
 */

template <typename CoeffType>
class polynomial;

template <typename CoeffType>
polynomial<CoeffType> operator * (polynomial<CoeffType> const& p1, polynomial<CoeffType> const& p2)
{
    polynomial<CoeffType> result;
    //            bool overflow = false;
    std::size_t max_order = p1.max_order;
    for(std::size_t je1 = 0; je1 < max_order; ++je1)
        for(std::size_t he1 = 0; he1 < max_order; ++he1)
        {
            for(std::size_t je2 = 0; je2 < max_order - je1; ++je2)
                for(std::size_t he2 = 0; he2 < max_order - he1; ++he2)
                    result.coeffs[ (je1+je2)*max_order + he1+he2 ] += p1.coeffs[je1*max_order+he1] * p2.coeffs[je2*max_order+he2];
            
            //                    // Overflow check
            //                    for(std::size_t je2 = max_order - je1; je2 < max_order; ++je2)
            //                        for(std::size_t he2 = max_order - he1; he2 < max_order; ++ he2)
            //                            if( coeffs[je2*max_order+he2] != CoeffType(0))
            //                                overflow = true;
        }
    //            if (overflow)
    //                std::cerr<<"WARNING: polynomial overflow -> truncated!"<<std::endl;
    
    return result;
}

template <typename CoeffType>
polynomial<CoeffType> operator * (polynomial<CoeffType> const& p, monomial<CoeffType> const& m)
{
    polynomial<CoeffType> r;
    std::size_t max_order = r.max_order;
    for(std::size_t je = 0; je < max_order-m.j_exp; ++je)
        for(std::size_t he = 0; he < max_order-m.h_exp; ++he)
            r(je+m.j_exp,he+m.h_exp) = p(je,he) * m.coeff;
    return r;
}
#endif //VLI__POLYNOMIAL
