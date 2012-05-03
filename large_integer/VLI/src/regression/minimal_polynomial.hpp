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
#ifdef _OPENMP
#include <omp.h>
#endif
namespace hp2c
{

/**
  * A monomial class
  */

template <typename CoeffType>
struct monomial
{
    std::size_t j_exp;
    std::size_t h_exp;
    CoeffType coeff;

    /**
      * Constructor: Creates a monomial 1*J^j_exp*h^h_exp
      */
    explicit monomial(unsigned int j_exp = 0, unsigned int h_exp = 0)
        : j_exp(j_exp), h_exp(h_exp), coeff(1)
    {
    }

    monomial& operator *= (CoeffType const& c)
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

/**
  * Multiplication with some factor of arbitrary type T
  * (for which the monomial class has to provide a *= operator)
  * e.g. T = int
  */
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

//
// *****************************************************************************
//

/**
  * A polynomial class
  */

template <typename CoeffType, unsigned int Order>
class polynomial;

/**
  * Multiplication of two polynomials
  */
template <typename CoeffType, unsigned int Order>
polynomial<CoeffType,2*Order> operator * (polynomial<CoeffType,Order> const& p1, polynomial<CoeffType,Order> const& p2)
{
    std::size_t max_order = p1.max_order;
    std::vector<CoeffType> result_coeffs(4*max_order*max_order);
    for(std::size_t je1 = 0; je1 < max_order; ++je1){
        for(std::size_t he1 = 0; he1 < max_order; ++he1){
            for(std::size_t je2 = 0; je2 < max_order; ++je2){
               for(std::size_t he2 = 0; he2 < max_order; ++he2){
                     result_coeffs[ 2*(je1+je2)*max_order + he1+he2 ] += p1.coeffs[je1*max_order+he1] * p2.coeffs[je2*max_order+he2];
               }
            }

        }
    }
    return polynomial<CoeffType,2*Order>(result_coeffs);
}

/**
  * Multiplication of a polynomial with a monomial
  */
template <typename CoeffType, unsigned int Order, typename T>
polynomial<CoeffType,Order> operator * (polynomial<CoeffType,Order> const& p, monomial<T> const& m)
{
    polynomial<CoeffType,Order> r;
    std::size_t max_order = r.max_order;
    for(std::size_t je = 0; je < max_order-m.j_exp; ++je)
        for(std::size_t he = 0; he < max_order-m.h_exp; ++he)
            r(je+m.j_exp,he+m.h_exp) = p(je,he)  * m.coeff;
    return r;
}

template <typename CoeffType, unsigned int Order>
class polynomial
{
    private:
        std::vector<CoeffType> coeffs;
    public:
        enum {max_order = Order};

        friend polynomial<CoeffType,2*Order> operator *<> (polynomial<CoeffType,Order> const&, polynomial<CoeffType,Order> const&);

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

        /**
          * Copy constructor
          */
        polynomial(polynomial const& p)
            :coeffs(p.coeffs)
        {
        }

        explicit polynomial(std::vector<CoeffType> const& v)
            : coeffs(v)
        {
            assert(v.size() == max_order*max_order);
        }

        /**
          * Assignment operator
          */
        polynomial& operator = (polynomial p)
        {
            swap(*this,p);
            return *this;
        }

        /**
          * Swap function
          */
        friend void swap(polynomial& p1, polynomial& p2)
        {
            swap(p1.coeffs,p2.coeffs);
        }

        /**
          * Prints polynomial to ostream o
          */
        void print(std::ostream& o) const
        {
            for(std::size_t je = 0; je < max_order; ++je)
                for(std::size_t he = 0; he < max_order; ++he)
                    if( coeffs[je*max_order+he] != CoeffType(0))
                        o<<"+"<<coeffs[je*max_order+he]<<"*J^"<<je<<"*h^"<<he<<std::endl;
        }

        /**
          * Plus assign with a polynomial
          */
        polynomial& operator += (polynomial const& p)
        {
            typename std::vector<CoeffType>::iterator it    = coeffs.begin();
            typename std::vector<CoeffType>::iterator end   = coeffs.end();
            typename std::vector<CoeffType>::const_iterator p_it  = p.coeffs.begin();
            typename std::vector<CoeffType>::const_iterator p_end = p.coeffs.end();
            while( it != end && p_it != p_end)
            {
                *it += *p_it;
                ++it;
                ++p_it;
            }
            return *this;
        }

        /**
          * Plus assign with a coefficient
          * Adds the coefficient to the 0th order element of the polynomial
          */
        polynomial& operator += (CoeffType c)
        {
            coeffs[0] += c;
            return *this;
        }

        /**
          * Plus assign with a monomial
          */
        template <typename T>
        polynomial& operator += (monomial<T> const& m)
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

        /**
          * Multiplies assign with coefficient
          * Mutliplies all elements the argument
          */
        template <typename T>
        polynomial& operator *= (T const& c)
        {
            for(typename std::vector<CoeffType>::iterator it = coeffs.begin(); it != coeffs.end(); ++it)
                *it *= c;
            return *this;
        }

        /**
         * Comparison with an CoeffType (0th order)
         */

        bool operator == (CoeffType const& c) const
        {
            if(coeffs[0] == c)
            {
                bool all_zero = true;
                for(typename std::vector<CoeffType>::const_iterator it = coeffs.begin(); it != coeffs.end(); ++it)
                    all_zero = all_zero && (*it == 0);
                return all_zero;
            }
            return false;
        }

        /**
         * Access coefficient of monomial J^j_exp*h^h_exp
         */
        inline CoeffType const& operator ()(unsigned int j_exp, unsigned int h_exp) const
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


/**
  * Stream operator
  */
template <typename CoeffType, unsigned int Order>
std::ostream& operator <<(std::ostream& o, polynomial<CoeffType,Order> const& p)
{
    p.print(o);
    return o;
}

/**
  * Multiplication of a monomial with a polynomial
  */
template <typename CoeffType, unsigned int Order, typename T>
polynomial<CoeffType,Order> operator * (monomial<T> const& m,polynomial<CoeffType,Order> const& p)
{
    return p * m;
}

/**
  * Multiplication of a polynomial with a factor
  */
template <typename CoeffType, unsigned int Order, typename T>
polynomial<CoeffType,Order> operator * (polynomial<CoeffType,Order> p, T const& c)
{
    p *= c;
    return p;
}

template <typename CoeffType, unsigned int Order, typename T>
polynomial<CoeffType,Order> operator * (T const& c, polynomial<CoeffType,Order> const& p)
{
    return p * c;
}

template <typename CoeffType, unsigned int Order>
polynomial<CoeffType,2*Order>  inner_product(std::vector<polynomial<CoeffType,Order> > const& a, std::vector<polynomial<CoeffType,Order> >const& b)
{
    assert( a.size() == b.size() );


#ifdef _OPENMP
    std::vector < polynomial<CoeffType,2*Order> > result(omp_get_max_threads());
    #pragma omp parallel for
    for(int i=0; i < static_cast<int>(a.size());++i){
        result[omp_get_thread_num()] += a[i]*b[i];
    }
    
    for(int i=1; i < omp_get_max_threads();++i){
        result[0] += result[i];
    }
    
    return result[0];
#else
    polynomial<CoeffType,2*Order> result;
    typename std::vector<polynomial<CoeffType,Order> >::const_iterator it(a.begin()),  it_b(b.begin());
    while( it != a.end() )
    {
        result += *it * *it_b;
        ++it;
        ++it_b;
    }

    return result;
#endif
}

}

#endif //HP2C__MINIMAL_POLYNOMIAL
