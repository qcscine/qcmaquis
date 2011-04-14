/*****************************************************************************
 *
 * A minimal implementation of a polynomial in two symbolic variables 'J','h'.
 *
 * (C) 2010 by Andreas Hehn (hehn@phys.ethz.ch)
 *
 *****************************************************************************/

#ifndef HP2C__MINIMAL_POLYNOMIAL
#define HP2C__MINIMAL_POLYNOMIAL

#include <vector>
#include <iostream>

namespace hp2c
{

/**
  * A monomial class
  */

template <typename CoeffType>
struct monomial
{
    unsigned int j_exp;
    unsigned int h_exp;
    CoeffType coeff;

    /**
      * Constructor: Creates a monomial 1*J^j_exp*h^h_exp
      */
    monomial(unsigned int j_exp = 0, unsigned int h_exp = 0)
        : j_exp(j_exp), h_exp(h_exp), coeff(1)
    {
    }

    monomial& operator *= (CoeffType const& c)
    {
        coeff *= c;
        return *this;
    }
};

template <typename CoeffType>
monomial<CoeffType> operator * (monomial<CoeffType> m, CoeffType const& c)
{
    m *= c;
    return m;
}

template <typename CoeffType>
monomial<CoeffType> operator * (CoeffType const& c, monomial<CoeffType> const& m)
{
    return m*c;
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

template <typename CoeffType>
class polynomial
{
    private:
        std::vector<CoeffType> coeffs;
    public:
        enum {max_order = POLYNOMIAL_MAX_ORDER};

        friend polynomial operator *<> (polynomial const&, polynomial const&);

        /**
          * Constructor: Creates a polynomial
          * p(J,h) = c_0 + c_{1,0} * J^1 * h^0 + c_{2,0} * J^2 * h^0 + ...
          *             + c_{max_order-1,max_order-1} * (J*h)^(max_order-1)
          * where all c_* are of type CoeffType and set to 0.
          */
        polynomial()
            : coeffs(max_order*max_order,CoeffType(0))
        {
        }

        polynomial(polynomial const& p)
            :coeffs(p.coeffs)
        {
        }

        polynomial& operator = (polynomial p)
        {
            swap(p);
            return *this;
        }

        void swap(polynomial& p)
        {
            coeffs.swap(p.coeffs);
        }

        void print(std::ostream& o) const
        {
            for(std::size_t je = 0; je < max_order; ++je)
                for(std::size_t he = 0; he < max_order; ++he)
                    if( coeffs[je*max_order+he] != CoeffType(0))
                        o<<"+"<<coeffs[je*max_order+he]<<"*J^"<<je<<"*h^"<<he;
        }

        polynomial& operator += (polynomial const& p)
        {
            typename std::vector<CoeffType>::iterator it    = coeffs.begin();
            typename std::vector<CoeffType>::iterator end   = coeffs.end();
            typename std::vector<CoeffType>::const_iterator p_it  = p.coeffs.begin();
            typename std::vector<CoeffType>::const_iterator p_end = p.coeffs.end();
            while( it != end, p_it != p_end)
            {
                *it += *p_it;
                ++it;
                ++p_it;
            }
            return *this;
        }

        polynomial& operator += (CoeffType const& c)
        {
            coeffs[0] += c;
            return *this;
        }

        polynomial& operator += (monomial<CoeffType> const& m)
        {
            assert(m.j_exp < max_order);
            assert(m.h_exp < max_order);
            coeffs[m.j_exp*max_order+m.h_exp] += m.coeff;
            return *this;
        }


//        polynomial& operator *= (monomial<CoeffType> const& m)
//        {
//            polynomial p;
//            for(std::size_t je = 0; je < max_order; ++je)
//                for(std::size_t he = 0; he < max_order; ++he)
//                    p.coeffs[ (je+m.j_exp)*max_order + he+m.h_exp ] = m.coeff * this->coeffs[ je*max_order + he ];
//            swap(p);
//            return *this;
//        }

        polynomial& operator *= (CoeffType const& c)
        {
            for(typename std::vector<CoeffType>::iterator it = coeffs.begin(); it != coeffs.end(); ++it)
                *it *= c;
            return *this;
        }

        /**
         * Access coefficient of monomial J^j_exp*h^h_exp
         */

        inline CoeffType operator ()(unsigned int j_exp, unsigned int h_exp) const
        {
            assert(j_exp < max_order);
            assert(h_exp < max_order);
            return coeffs[j_exp*max_order+h_exp];
        }
        
        /**
         * Access coefficient of monomial J^j_exp*h^h_exp
         */

        inline CoeffType& operator ()(unsigned int j_exp, unsigned int h_exp)
        {
            assert(j_exp < max_order);
            assert(h_exp < max_order);
            return coeffs[j_exp*max_order+h_exp];
        }
};

template <typename CoeffType>
std::ostream& operator <<(std::ostream& o, polynomial<CoeffType> const& p)
{
    p.print(o);
    return o;
}

template <typename CoeffType>
polynomial<CoeffType> operator * (monomial<CoeffType> const& m,polynomial<CoeffType> const& p)
{
    return p * m;
}

template <typename CoeffType>
polynomial<CoeffType> operator * (polynomial<CoeffType> p, CoeffType const& c)
{
    p *= c;
    return p;
}

template <typename CoeffType>
polynomial<CoeffType> operator * (CoeffType const& c, polynomial<CoeffType> const& p)
{
    return p * c;
}

}

#endif //HP2C__MINIMAL_POLYNOMIAL
