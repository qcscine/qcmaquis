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
#include "boost/swap.hpp"
#include "monome/polynome_cpu.h"
#include "monome/polynome_gpu.h"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "function_hooks/vli_vector_polynomial_gpu_function_hooks.hpp"

#include <ostream>


namespace vli
{



    template<class polynomial_cpu> 
    class vector_polynomial_cpu : public std::vector<polynomial_cpu>{ 
    private:
        typedef typename polynomial_cpu::value_type vli_value_type; 
        typedef typename std::size_t size_t;
        enum {max_order_poly = polynomial_cpu::max_order };
        enum {vli_size   = polynomial_cpu::size }; 
        enum {OffSet = max_order_poly*max_order_poly*vli_size}; 
    public:
    
    vector_polynomial_cpu(size_t size = 1)
    :std::vector<polynomial_cpu>(size)
    {
    }
    //copy and assignemant are done by the std vector    
    };
    template <class BaseInt, int Size, int Order>
    vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> > 
    inner_product( vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v1, 
                   vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v2){
        assert(v1.size() == v2.size());
        size_t size_v = v1.size();
        vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  res;
        
        for(size_t i=0 ; i < size_v ; i++){
            res[0] += v1[i]*v2[i]; 
        }

        return res;
    }

    template<class BaseInt, int Size, int Order > // the << cpu and gpu should be commun
 	std::ostream & operator<<(std::ostream & os, vector_polynomial_cpu< polynomial_cpu< vli_cpu<BaseInt, Size>, Order > >   & v)
    {
        os << "---------- CPU ----------" << std::endl;
        
        for(std::size_t i = 0; i < v.size(); i++)
            os << v[i] << std::endl;
        return os;
    }
}

#endif //VLI_VECTOR_POLYNOME_GPU_H
