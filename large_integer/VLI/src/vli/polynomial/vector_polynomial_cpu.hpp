/*
 *  vector_polynomial_cpu.h
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 10.07.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */

#ifndef VLI_VECTOR_POLYNOME_CPU_H
#define VLI_VECTOR_POLYNOME_CPU_H
#include "vli/polynomial/polynomial_cpu.hpp"
#include "vli/polynomial/polynomial_gpu.hpp"
#include "vli/polynomial/vector_polynomial_gpu.hpp"
#include "vli/vli_cpu.hpp"
#include "vli/detail/vli_size_param.hpp"
#include <vector>
#include <ostream>
#include <cassert>
#include <boost/swap.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif //_OPENMP

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
    
    
    template <class BaseInt, std::size_t Size, unsigned int Order>
    polynomial_cpu<vli_cpu<BaseInt, Size>, Order> 
    inner_product( vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v1, 
                   vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v2){
        assert(v1.size() == v2.size());
        std::size_t size_v = v1.size();
        
#ifdef _OPENMP
        polynomial_cpu<vli_cpu<BaseInt, Size>, Order>  res[omp_get_max_threads()];
        polynomial_cpu<vli_cpu<BaseInt, Size>, Order>  res_gpu_to_cpu;
        
        /**
        * 
        * Quick Mixmode
        */
        std::size_t split = static_cast<std::size_t>(v1.size()*0.5);
        polynomial_gpu<vli_gpu<BaseInt, Size>, Order>  resGpu;

        vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> > v1Gpu,v2Gpu;
        v1Gpu.resize(split);
        v2Gpu.resize(split);
        v1Gpu.copy_vec_vli_to_gpu(v1);
        v2Gpu.copy_vec_vli_to_gpu(v2);
        resGpu = inner_product(v1Gpu,v2Gpu);

        #pragma omp parallel for
        for(std::size_t i=split+1 ; i < size_v ; ++i){
            res[omp_get_thread_num()] += v1[i]*v2[i];
        }

        for(int i=1; i < omp_get_max_threads(); ++i)
            res[0]+=res[i];
        
        return res[0];
#else //_OPENMP
        polynomial_cpu<vli_cpu<BaseInt, Size>, Order>  res;
        for(std::size_t i=0 ; i < size_v ; ++i)
            res += v1[i]*v2[i];

        res += polynomial_cpu<vli_cpu<vli_value_type, Vli::size>, Order >(resGpu);


        return res;
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
