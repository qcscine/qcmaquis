#ifndef HP2C__SIMPLE_POLYNOMIAL_HPP
#define HP2C__SIMPLE_POLYNOMIAL_HPP
#include <vector>
#include <stdexcept>
#include <iostream>

namespace series_expansion
{


template <typename T>
class simple_polynomial
{
    private:
        std::vector<T> orders;
    public:
        typedef T value_type;
        typedef typename std::vector<T>::iterator       iterator;
        typedef typename std::vector<T>::const_iterator const_iterator;


        simple_polynomial()
            :orders(0,T(0))
            {}
        explicit simple_polynomial(T const& t)			// This constructor sets the 0th order value
            :orders(1,T(t))
            {}
        ~simple_polynomial()
            {}

        unsigned int size() const
        {
            return orders.size();
        }

        void truncate(unsigned int order)
        {
            if(order+1 < orders.size())
                orders.resize(order+1);
        }

		const T operator [] (unsigned int i) const
        {
            if(i< orders.size())
                return orders[i];
            else
                return T(0);
        }
        
		T& operator [] (unsigned int i)
        {
            if(i >= orders.size())
                orders.resize(i+1,T(0));
            return orders[i];
        }
		
		simple_polynomial& operator += (simple_polynomial const& rhs)
        {
            if(rhs.size() > this->size())
                orders.resize(rhs.size(),T(0));
            for(unsigned int i=0; i<rhs.size(); ++i)
                orders[i] += rhs.orders[i];
            return *this;
        }

        simple_polynomial& operator -= (simple_polynomial const& rhs)
        {
            if(rhs.size() > this->size())
                orders.resize(rhs.size(),T(0));
            for(unsigned int i=0; i<rhs.size(); ++i)
                orders[i] -= rhs.orders[i];
            return *this;
        }

		simple_polynomial& operator *= (T const& rhs)
        {
            for(iterator it(orders.begin()); it < orders.end(); ++it)
                *it *= rhs;
            return *this;
        }

        simple_polynomial& operator *= (simple_polynomial const& rhs)
        {
            std::runtime_error("Not implemented yet.");
            return *this;
        }

        simple_polynomial& operator /= (T const& rhs)
        {
            for(iterator it(orders.begin()); it < orders.end(); ++it)
                *it /= rhs;
            return *this;
        }
};

template <typename T>
inline const simple_polynomial<T> operator + (simple_polynomial<T> a, simple_polynomial<T> const& b)
{
    a+=b;
    return a;
}

template <typename T>
inline const simple_polynomial<T> operator - (simple_polynomial<T> a, simple_polynomial<T> const& b)
{
    a-=b;
    return a;
}

template <typename T>
inline const simple_polynomial<T> operator * (T const& t, simple_polynomial<T> a)
{
    a *= t;
    return a;
}

template <typename T>
inline const simple_polynomial<T> operator * (simple_polynomial<T> a, T const& t)
{
    a *= t;
    return a;
}

template <typename T>
inline const simple_polynomial<T> operator * (simple_polynomial<T> a, simple_polynomial<T> const& b)
{
    a*=b;
    return a;
}

template <typename T>
inline const simple_polynomial<T> operator / (simple_polynomial<T> a, T const& t)
{
    a /= t;
    return a;
}

template <typename T>
std::ostream& operator << (std::ostream &o, simple_polynomial<T> const& a)
{
    for(unsigned int i=0; i< a.size(); ++i)
    {
        if(a[i] < 0)
            o<<a[i]<<"*l^"<<i;
        else
            o<<"+"<<a[i]<<"*l^"<<i;
    }
    return o;
}

}

#endif //HP2C__SIMPLE_POLYNOMIAL_HPP
