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
#include "vli/polynomial/polynomial_cpu.hpp"
#include "vli/vli_cpu.hpp"
#include "vli/detail/kernels_cpu_gpu.hpp"
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

    template<class polynomial_cpu>
    class vector_polynomial_cpu : public std::vector<polynomial_cpu>{ 
    private:
        typedef typename polynomial_cpu::value_type::value_type vli_value_type;
        enum {max_order_poly = polynomial_cpu::max_order };
        enum {vli_size   = polynomial_cpu::value_type::size }; 
        enum {element_offset = max_order_poly*max_order_poly*vli_size}; 
    public:
    
    vector_polynomial_cpu(std::size_t size = 1)
    :std::vector<polynomial_cpu>(size)
    {
    }
    //copy and assignemant are done by the std vector   
            
    }; //end class
    
namespace detail
{    
#ifdef _OPENMP
    template <class BaseInt, std::size_t Size, unsigned int Order>
    polynomial_cpu<vli_cpu<BaseInt, Size>, Order> 
    inner_product_openmp( vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v1, 
                   vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v2){
        assert(v1.size() == v2.size());
        std::size_t size_v = v1.size();
    
        /**
        * the first version used an array:
        * polynomial_cpu<vli_cpu<BaseInt, Size>, Order>  res[omp_get_max_threads()];
        * GCC and C99 allow an array's size to be determined at run time. This extension is not permitted in standard C++.
        * so vector !
        */

        std::vector<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> > res(omp_get_max_threads()); 
       
        #pragma omp parallel for
        for(std::size_t i=0 ; i < size_v ; ++i){
            res[omp_get_thread_num()] += v1[i]*v2[i];
        }

        for(int i=1; i < omp_get_max_threads(); ++i)
            res[0]+=res[i];
        return res[0];
    }
#endif
    
    template <class BaseInt, std::size_t Size, unsigned int Order>
    polynomial_cpu<vli_cpu<BaseInt, Size>, Order> 
    inner_product_plain( vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v1, 
                   vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v2){
        assert(v1.size() == v2.size());
        std::size_t size_v = v1.size();
        polynomial_cpu<vli_cpu<BaseInt, Size>, Order>  res;
        for(std::size_t i=0 ; i < size_v ; ++i)
            res += v1[i]*v2[i];
        return res;
    }

    template <class BaseInt, std::size_t Size, unsigned int Order>
    polynomial_cpu<vli_cpu<BaseInt, Size>, Order> 
    inner_product_accp( vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  & v1, 
                   vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  & v2){
        assert(v1.size() == v2.size());
        typedef typename polynomial_cpu<vli_cpu<BaseInt,Size>,Order>::exponent_type exponent_type;
        std::size_t size_v = v1.size();
        vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >   vres(size_v);
        polynomial_cpu<vli_cpu<BaseInt, Size>, Order> res, res0;
        vli_cpu<BaseInt, Size> resvli; 
        bool result_is_negative;
        int pencil1, pencil2, pencil3;
        BaseInt r[2] = {0,0};	//for local block calculation


        // want 
        // we do the next following line, I tried to unroll/inline everything, but it came slower !

        // C - Tim you're stupid ! RTFM !

        #pragma omp acc_region
        #pragma omp acc_loop
        for(std::size_t i=0 ; i < size_v ; ++i)
        {
            for(exponent_type je1 = 0; je1 < Order; ++je1)
            {
                for(exponent_type je2 = 0; je2 < Order - je1; ++je2)
                {
                    for(exponent_type he1 = 0; he1 < Order; ++he1)
                    {
                        for(exponent_type he2 = 0; he2 < Order - he1; ++he2)
                        {  
                            // we want
                            // res.coeffs_[(je1+je2)*Order + he1+he2 ] += v1[i].coeffs_[je1*Order+he1] * v2[i].coeffs_[je2*Order+he2];
                            // we do resvli = v1[i].coeffs_[je1*Order+he1] * v2[i].coeffs_[je2*Order+he2]
                            // and res.coeffs_[(je1+je2)*Order + he1+he2 ] += resvli
                            // OK I can not do more for the cray compiler
                            pencil1=je1*Order+he1;
                            pencil2=je2*Order+he2;
                            pencil3=(je1+je2)*Order + he1+he2;

                            result_is_negative = static_cast<bool>((v1[i].coeffs_[pencil1].data_[Size-1] ^ v2[i].coeffs_[pencil2].data_[Size-1]) >> data_bits<BaseInt>::value);
                          
                            if(result_is_negative)// test if, for the negative case ...
                            {
                                v1[i].coeffs_[pencil1].negate(); // - to +
                                kernels_multiplication_classic_truncate<BaseInt,Size>(&resvli[0],&v1[i].coeffs_[pencil1].data_[0], &v2[i].coeffs_[pencil2].data_[0]);
                                v1[i].coeffs_[pencil1].negate(); // + to -
                            }else{                          
                              kernels_multiplication_classic_truncate<BaseInt,Size>(&resvli[0],&v1[i].coeffs_[pencil1].data_[0], &v2[i].coeffs_[pencil2].data_[0]);                            
                            }
    //                        res.coeffs_[(je1+je2)*Order + he1+he2 ] += resvli;
                         
                            for( std::size_t i = 0; i < Size-1 ; ++i){
                                res.coeffs_[pencil3].data_[i]   += resvli[i]; 
                                res.coeffs_[pencil3].data_[i+1] += res.coeffs_[pencil3].data_[i]  >> data_bits<BaseInt>::value; //carry bit
                                res.coeffs_[pencil3].data_[i]   &= data_mask<BaseInt>::value; // Remove the carry bit
                            }
     
                      	    res.coeffs_[pencil3].data_[Size-1 ] += resvli[Size-1];
                            res.coeffs_[pencil3].data_[Size-1 ] &= base<BaseInt>::value + data_mask<BaseInt>::value;

                            resvli[0]=0;
                            resvli[1]=0;
                            resvli[2]=0;
                        }
                    }
                }    
            }
        }

        return res;
    }
} // end namespace detail

    template <class BaseInt, std::size_t Size, unsigned int Order>
    inline polynomial_cpu<vli_cpu<BaseInt, Size>, Order> 
    inner_product( vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v1, 
                   vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v2){
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

    template<class BaseInt, std::size_t Size, unsigned int Order > // the << cpu and gpu should be commun
 	std::ostream & operator<<(std::ostream & os, vector_polynomial_cpu< polynomial_cpu< vli_cpu<BaseInt, Size>, Order > > const& v)
    {        
        for(std::size_t i = 0; i < v.size(); ++i)
            os << v[i] << std::endl;
        return os;
    }
}

#endif //VLI_VECTOR_POLYNOME_GPU_H
