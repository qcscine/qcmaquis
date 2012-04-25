/*
 *  vector_polynomial_cpu.h
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 10.07.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */

#ifndef VLI_VECTOR_POLYNOME_CPU_H
#define VLI_VECTOR_POLYNOME_CPU_H
#include "vli/vli_config.h"
#include "vli/polynomial/polynomial.hpp"
#include "vli/vli_cpu.h"
#include <vector>
#include <ostream>
#include <cassert>
#include <boost/swap.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif //_OPENMP

#ifdef VLI_COMPILE_GPU
#include "vli/detail/inner_product_gpu_booster.hpp"
#endif //VLI_COMPILE_GPU

namespace vli
{

template <class BaseInt, std::size_t Size>
class vli_cpu;

template<class Polynomial>
class vector_polynomial : public std::vector<Polynomial>{ 
  private:
//    typedef typename Polynomial::value_type::value_type vli_value_type;
    
//    static typename exponent_type<Polynomial>::type max_order_poly;
//    enum {vli_size   = polynomial_cpu::value_type::size }; 
    static std::size_t const element_offset = Polynomial::max_order*Polynomial::max_order*Polynomial::value_type::size;

  public:
    vector_polynomial(std::size_t size = 1)
    :std::vector<Polynomial>(size) {
    }
    //copy and assignemant are done by the std vector
}; //end class

namespace detail
{    
#ifdef _OPENMP
template <class BaseInt, std::size_t Size, unsigned int Order>
polynomial<vli_cpu<BaseInt, 2*Size>, 2*Order>
inner_product_openmp( vector_polynomial<polynomial<vli_cpu<BaseInt, Size>, Order> >  const& v1, 
                      vector_polynomial<polynomial<vli_cpu<BaseInt, Size>, Order> >  const& v2){
    assert(v1.size() == v2.size());
    std::size_t size_v = v1.size();

    std::vector<polynomial<vli_cpu<BaseInt, 2*Size>, 2*Order> > res(omp_get_max_threads()); 
   
    #pragma omp parallel for
    for(long i=0 ; i < size_v ; ++i){
        res[omp_get_thread_num()] += v1[i]*v2[i];
    }

    for(int i=1; i < omp_get_max_threads(); ++i)
        res[0]+=res[i];

    return res[0];

}
#endif

template <class BaseInt, std::size_t Size, unsigned int Order>
polynomial<vli_cpu<BaseInt, 2*Size>, 2*Order> 
inner_product_plain( vector_polynomial<polynomial<vli_cpu<BaseInt, Size>, Order> >  const& v1, 
               vector_polynomial<polynomial<vli_cpu<BaseInt, Size>, Order> >  const& v2){
    assert(v1.size() == v2.size());
    std::size_t size_v = v1.size();
    polynomial<vli_cpu<BaseInt, 2*Size>, 2*Order>  res;
    for(std::size_t i=0 ; i < size_v ; ++i)
        res += v1[i]*v2[i];
    return res;
}
} // end namespace detail


template <class VectorPolynomial>
struct inner_product_result_type {
};

template <class Coeff, unsigned int Order>
struct inner_product_result_type< vector_polynomial<polynomial<Coeff,Order> > > {
    typedef typename polynomial_multiply_result_type<polynomial<Coeff,Order> >::type type;
};


template <class Coeff, unsigned int Order>
inline typename inner_product_result_type<vector_polynomial<polynomial<Coeff,Order> > >::type
inner_product(
          vector_polynomial<polynomial<Coeff, Order> > const& v1
        , vector_polynomial<polynomial<Coeff, Order> > const& v2
        ) {
    // TODO this is a little dirty and could be done better
#ifdef _OPENMP
#ifdef VLI_USE_GPU
    return detail::inner_product_openmp_gpu(v1,v2);
#else //VLI_USE_GPU
    return detail::inner_product_openmp(v1,v2);
#endif //VLI_USE_GPU
#else //_OPENMP
#ifdef VLI_USE_GPU
    return detail::inner_product_gpu(v1,v2);
#else //VLI_USE_GPU
    return detail::inner_product_plain(v1,v2);
#endif //VLI_USE_GPU
#endif //_OPENMP
}

template<class Polynomial, unsigned int Order >
std::ostream & operator<<(std::ostream & os, vector_polynomial<Polynomial> const& v)
{        
    for(std::size_t i = 0; i < v.size(); ++i)
        os << v[i] << std::endl;
    return os;
}
}

#endif //VLI_VECTOR_POLYNOME_GPU_H
