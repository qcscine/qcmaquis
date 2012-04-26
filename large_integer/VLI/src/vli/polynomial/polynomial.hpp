#ifndef VLI_POLYNOME_CPU_H
#define VLI_POLYNOME_CPU_H
#include "vli/function_hooks/vli_polynomial_cpu_function_hooks.hpp"
#include "vli/polynomial/monomial.hpp"
#include "vli/vli_cpu.h"

#include <boost/swap.hpp>
#include <boost/type_traits/is_fundamental.hpp>
#include <ostream>
#include <cassert>

namespace vli
{	    

template <class Coeff, unsigned int Order>
class polynomial;


//------------------------------------------------------------------------ 
// Polynomial traits
//------------------------------------------------------------------------

template <typename Polynomial>
struct polynomial_multiply_result_type {
};

template <typename Coeff, unsigned int Order>
struct polynomial_multiply_result_type<polynomial<Coeff,Order> > {
    typedef polynomial<Coeff,2*Order> type;
};

template <typename BaseInt, std::size_t Size, unsigned int Order>
struct polynomial_multiply_result_type<polynomial<vli_cpu<BaseInt,Size>,Order> > {
    typedef polynomial<vli_cpu<BaseInt,2*Size>,2*Order> type;
};

template <typename Polynomial>
struct exponent_type{
    typedef typename Polynomial::exponent_type type;
};

template <typename T>
bool is_zero(T t) {
    return t == 0;
}

template <typename T>
void negate_inplace(T& t) {
    t = -t;
}

template <typename T, typename T2, typename T3>
void muladd(T& t, T2 const& t2, T3 const& t3) {
    t += t2 * t3;
}


namespace detail {

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

    template <class Coeff, unsigned int Order>
    void init(polynomial<Coeff,Order>& p, boost::false_type dummy) {
        // This is a non fundamental type -> it will be default constructed
    }
    
    template <class Coeff, unsigned int Order>
    void init(polynomial<Coeff,Order>& p, boost::true_type dummy) {
        // This is a fundamental type (e.g. unsigned int, double,...) -> we have to initalize
        for(typename polynomial<Coeff,Order>::exponent_type i=0; i<Order*Order;++i)
            p.coeffs_[i]=Coeff();
    }

    template <class Operation, class Coeff, unsigned int Order>
    void additive_op_assign(polynomial<Coeff,Order>& p, polynomial<Coeff,Order> const& p2, Operation op) {
        for(typename polynomial<Coeff,Order>::exponent_type i=0; i<Order*Order; ++i)
            op(p.coeffs_[i],p2.coeffs_[i]);
    }

    template <class Operation, class Coeff, unsigned int Order>
    void additive_op_assign(polynomial<Coeff,Order>& p, monomial<Coeff> const& m, Operation op) {
        op(p(m.j_exp_,m.h_exp_), m.coeff_);
    }
    
    template <class Operation, class Coeff, unsigned int Order>
    void additive_op_assign(polynomial<Coeff,Order>& p, Coeff const& c, Operation op) {
        op(p(0,0),c);
    }

    template <class Operation, class Coeff, unsigned int Order>
    void additive_op_assign(polynomial<Coeff,Order>& p, int a, Operation op) {
        op(p(0,0),a);
    }
    
    template <class Operation, class Coeff, unsigned int Order>
    void multiplicative_op_assign(polynomial<Coeff,Order>& p, int a, Operation op) {
        for(typename polynomial<Coeff,Order>::exponent_type i=0; i< Order*Order; ++i)
            op(p.coeffs_[i], a);
    }
    
    template <class Operation, class Coeff, unsigned int Order>
    void multiplicative_op_assign(polynomial<Coeff,Order>& p, Coeff const& c, Operation op) {
        for(typename polynomial<Coeff,Order>::exponent_type i=0; i< Order*Order; ++i)
            op(p.coeffs_[i], c);
    }
    
    template <class Coeff, unsigned int Order, class T>
    void multiply_assign_monomial(polynomial<Coeff,Order>& p, monomial<T> const& m) {
        typedef typename polynomial<Coeff,Order>::value_type value_type;
        for(std::ptrdiff_t je=static_cast<std::ptrdiff_t>(Order)-1-m.j_exp_; je >= 0; --je)
        {
            for(std::ptrdiff_t he=static_cast<std::ptrdiff_t>(Order)-1-m.h_exp_; he >= 0; --he)
                p(je+m.j_exp_, he+m.h_exp_) = p(je,he)*m.coeff_;
            for(std::ptrdiff_t he=static_cast<std::ptrdiff_t>(m.h_exp_)-1; he >=0; --he)
                p(je+m.j_exp_,he) = value_type(0);
        }

        for(std::ptrdiff_t je=static_cast<std::ptrdiff_t>(m.j_exp_)-1; je >=0; --je)
            for(std::ptrdiff_t he=static_cast<std::ptrdiff_t>(Order)-1; he >=0; --he)
                p(je,he) = value_type(0);
    }

    template <class Coeff, unsigned int Order>
    void negate(polynomial<Coeff,Order>& p) {
        for(typename polynomial<Coeff,Order>::exponent_type i = 0; i < Order*Order; ++i)
            negate_inplace(p.coeffs_[i]);
    }
    
    template <class Coeff, unsigned int Order>
    bool is_zero_helper(polynomial<Coeff,Order> const& p)
    {
        for(typename polynomial<Coeff,Order>::exponent_type i=0; i<Order*Order; ++i)
            if (!is_zero(p.coeffs_[i]))
                return false;
        return true;
    }

    template <class Polynomial>
    struct equal_helper {
        bool operator()(Polynomial const& p, Polynomial const& p2) {
            typedef typename Polynomial::exponent_type exponent_type;
            for(exponent_type i = 0; i < Polynomial::max_order*Polynomial::max_order; ++i)
                if (!(p.coeffs_[i] == p2.coeffs_[i]))
                    return false;
            return true;
        }
    };
    
    template <class Polynomial>
    struct polynomial_multiply_helper {
    };

    template <class Coeff, unsigned int Order>
    struct polynomial_multiply_helper<polynomial<Coeff,Order> > {
        typename polynomial_multiply_result_type<polynomial<Coeff,Order> >::type operator()(polynomial<Coeff,Order> const& p1, polynomial<Coeff,Order> const& p2) {
            typedef typename polynomial<Coeff,Order>::exponent_type exponent_type;
            typename polynomial_multiply_result_type<polynomial<Coeff,Order> >::type result;
            for(exponent_type je1 = 0; je1 < Order; ++je1)
                for(exponent_type je2 = 0; je2 < Order; ++je2)
                    for(exponent_type he1 = 0; he1 < Order; ++he1)
                        for(exponent_type he2 = 0; he2 < Order; ++he2)
                            muladd(result.coeffs_[(je1+je2)*2*Order + he1+he2 ], p1.coeffs_[je1*Order+he1],p2.coeffs_[je2*Order+he2]);
            return result;
        }
    };

} // end namespace detail

//------------------------------------------------------------------------ 
// Polynomial template class
//------------------------------------------------------------------------ 

template <class Coeff, unsigned int Order>          
void swap(polynomial<Coeff,Order>& p1, polynomial<Coeff,Order>& p2){
    using boost::swap;
    swap(p1.coeffs_,p2.coeffs_);
}

template <class Coeff, unsigned int Order>
class polynomial{
public:
    typedef unsigned int exponent_type;      // Type of the exponents (has to be the same type as Vli::size_type)
    typedef Coeff      value_type;
    typedef value_type coeff_type;
    static exponent_type const max_order = Order;
        
    polynomial() {
        detail::init(*this, typename boost::is_fundamental<Coeff>::type());
    }
    explicit polynomial(int i) {
        detail::init(*this, typename boost::is_fundamental<Coeff>::type());
        coeffs_[0] = value_type(i);
    }
    polynomial(polynomial const& p) { 
        for(exponent_type i=0; i<Order*Order;++i)
            coeffs_[i]=p.coeffs_[i];
    }
    template <class Coeff2>
    explicit polynomial(polynomial<Coeff2,Order> const& p) {
        for(exponent_type i=0; i<Order*Order; ++i)
            coeffs_[i] = static_cast<Coeff>(p.coeffs_[i]);
    }
    polynomial& operator = (polynomial p) {
        swap(*this,p);
        return *this;
    }
    
    polynomial& operator += (polynomial const& p) { 
        detail::additive_op_assign(*this, p, detail::plus_assign());
        return *this;
    }
    polynomial& operator += (monomial<Coeff> const& m) {
        detail::additive_op_assign(*this, m, detail::plus_assign());
        return *this;
    }
    polynomial& operator += (int a) {
        detail::additive_op_assign(*this, a, detail::plus_assign());
        return *this;
    }


    polynomial& operator -= (polynomial const& p) {
        detail::additive_op_assign(*this, p, detail::minus_assign());
        return *this;
    }
    polynomial& operator -= (monomial<Coeff> const& m) {
        detail::additive_op_assign(*this, m, detail::minus_assign());
        return *this;
    }
    polynomial& operator -= (int a) {
        detail::additive_op_assign(*this, a, detail::minus_assign());
        return *this;
    }
   

    template <class T>
    polynomial& operator *= (monomial<T> const& m) {
        detail::multiply_assign_monomial(*this,m);
        return *this;
    }
    polynomial& operator *= (Coeff const& c) {
        detail::multiplicative_op_assign(*this, c, detail::multiply_assign());
        return *this;
    }
    polynomial& operator *= (int a) {
        detail::multiplicative_op_assign(*this, a, detail::multiply_assign());
        return *this;
    }

    polynomial& operator /= (Coeff const& c) {
        detail::multiplicative_op_assign(*this, c, detail::devide_assign());
        return *this;
    }
    polynomial& operator /= (int a) {
        detail::multiplicative_op_assign(*this, a, detail::devide_assign());
        return *this;
    }

    polynomial operator - () const {
        polynomial r(*this);
        r.negate();
        return r;
    }
    
    void negate() {
        detail::negate(*this);
    }

    bool operator==(polynomial const& p) const {
        return detail::equal_helper<polynomial>()(*this,p);
    }
    bool is_zero() const {
        return detail::is_zero_helper(*this);
    }
    
    friend void swap<>(polynomial<Coeff,Order>& p1, polynomial<Coeff,Order>& p2);

    inline Coeff const& operator ()(exponent_type j_exp, exponent_type h_exp) const {
        assert(j_exp < max_order);
        assert(h_exp < max_order);
        return coeffs_[j_exp*max_order+h_exp];
    }
    
    inline Coeff& operator ()(exponent_type j_exp, exponent_type h_exp) {
        assert(j_exp < max_order);
        assert(h_exp < max_order);
        return coeffs_[j_exp*max_order+h_exp];
    }
    
    void print(std::ostream& os) const {
        for(exponent_type i = 0; i < Order ; i++) {
            for(exponent_type j = 0; j < Order ; j++) {
                if( coeffs_[i*Order+j] != Coeff(0)) {
                    if( !(coeffs_[i*Order+j] < Coeff(0)) )
                        os << "+";
                    os << coeffs_[i*Order+j]<<"*J^"<<i<<"*h^"<<j<<std::endl;
                }
            }
        }
    }

    coeff_type coeffs_[Order*Order];
};

template <class Coeff, unsigned int Order>
bool is_zero(polynomial<Coeff,Order> const& p) {
    return p.is_zero();
}

template <class Coeff, unsigned int Order>
void negate_inplace(polynomial<Coeff, Order>& p) {
    p.negate();
}

template <class Coeff, unsigned int Order>
polynomial<Coeff, Order> operator + (polynomial<Coeff,Order> a, polynomial<Coeff,Order> const& b) {
    a += b;
    return a;
}

template <class T, class Coeff, unsigned int Order>
polynomial<Coeff, Order> operator * (polynomial<Coeff, Order> p, monomial<T> const& m) {
    p *= m;
    return p;
}

template <class T, class Coeff, unsigned int Order>
polynomial<Coeff, Order> operator * (monomial<T> const& m,polynomial<Coeff, Order> const& p) {
    return p * m;
}

template <class Coeff, unsigned int Order>
polynomial<Coeff, Order> operator * (polynomial<Coeff, Order> p, Coeff const& c) {
    p *= c;
    return p;
}

template <class Coeff, unsigned int Order>
polynomial<Coeff, Order> operator * (polynomial<Coeff, Order> p, int c) {
    p *= c;
    return p;
}

template <class Coeff, unsigned int Order>
polynomial<Coeff, Order> operator * (Coeff const& c, polynomial<Coeff, Order> const& p) {
    return p * c;
}

template <class Coeff, unsigned int Order>
polynomial<Coeff, Order> operator * (int c, polynomial<Coeff, Order> const& p) {
    return p * c;
}

template <class Coeff, unsigned int Order>
polynomial<Coeff, Order> operator / (polynomial<Coeff,Order> p, int c) {
    p /= c;
    return p;
}

template <class Coeff, unsigned int Order>
polynomial<Coeff, Order> operator / (polynomial<Coeff,Order> p, Coeff const& c) {
    p /= c;
    return p;
}


template <class Coeff, unsigned int Order>
typename polynomial_multiply_result_type<polynomial<Coeff,Order> >::type operator * (polynomial<Coeff,Order> const& p1, polynomial<Coeff,Order> const& p2) {
    return detail::polynomial_multiply_helper<polynomial<Coeff,Order> >()(p1,p2);
}

template <class Coeff, unsigned int Order> 
std::ostream& operator<<(std::ostream& os, polynomial<Coeff, Order> const& p) {
    p.print(os);
    return os;
}


//------------------------------------------------------------------------ 
// Specializations for vli_cpu
//------------------------------------------------------------------------ 

template <class BaseInt, std::size_t Size>
class vli_cpu;

namespace detail {
    template <class BaseInt, std::size_t Size, unsigned int Order>
    struct equal_helper<polynomial<vli_cpu<BaseInt,Size>,Order> > {
        bool operator()(polynomial<vli_cpu<BaseInt,Size>,Order> const& p, polynomial<vli_cpu<BaseInt,Size>,Order> const& p2) {
            int n = memcmp((void*)&p.coeffs_[0],(void*)&p2.coeffs_[0],Order*Order*vli_cpu<BaseInt,Size>::size*sizeof(typename vli_cpu<BaseInt,Size>::value_type));
            return (0 == n);
        }
    };
}

} //end namespace

#endif
