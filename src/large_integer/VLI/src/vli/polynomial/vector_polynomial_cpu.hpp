/*
 *  vector_polynomial_cpu.h
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 10.07.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */

#ifndef VLI_VECTOR_POLYNOME_CPU_H
#define VLI_VECTOR_POLYNOME_CPU_H
#include <vector>
#include <omp.h>
#include <ostream>
#include <cassert>
#include "boost/swap.hpp"
#include "vli/polynomial/polynomial_cpu.hpp"
#include "vli/vli_cpu.hpp"
#include "vli/detail/vli_size_param.hpp"

namespace vli
{

    template <class BaseInt, int Size>
    class vli_cpu;

    template<class polynomial_cpu> 
    class vector_polynomial_cpu : public std::vector<polynomial_cpu>{ 
    private:
        typedef typename polynomial_cpu::value_type::value_type vli_value_type; 
        typedef typename std::size_t size_t;
        enum {max_order_poly = polynomial_cpu::max_order };
        enum {vli_size   = polynomial_cpu::value_type::size }; 
        enum {OffSet = max_order_poly*max_order_poly*vli_size}; 
    public:
    
    vector_polynomial_cpu(size_t size = 1)
    :std::vector<polynomial_cpu>(size)
    {
    }
    //copy and assignemant are done by the std vector   
            
    }; //end class
    
    
    template <class BaseInt, int Size, int Order>
    polynomial_cpu<vli_cpu<BaseInt, Size>, Order> 
    inner_product( vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v1, 
                   vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v2){
        assert(v1.size() == v2.size());
        size_t size_v = v1.size();
        vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> > inter;
        inter.resize(size_v);
        
        polynomial_cpu<vli_cpu<BaseInt, Size>, Order>  res[omp_get_max_threads()];
   
        /**
        * reduction is impossible on class type  
        * so type
        */
        #pragma omp parallel for
        for(std::size_t i=0 ; i < size_v ; ++i){            
            res[omp_get_thread_num()] += v1[i]*v2[i]; 
        }
        
        for(int i=1; i < omp_get_max_threads(); ++i)
            res[0]+=res[i];
        
        return res[0];
    }

    template<class BaseInt, int Size, int Order > // the << cpu and gpu should be commun
 	std::ostream & operator<<(std::ostream & os, vector_polynomial_cpu< polynomial_cpu< vli_cpu<BaseInt, Size>, Order > >   const& v)
    {        
        for(std::size_t i = 0; i < v.size(); i++)
            os << v[i] << std::endl;
        return os;
    }
}

#endif //VLI_VECTOR_POLYNOME_GPU_H
