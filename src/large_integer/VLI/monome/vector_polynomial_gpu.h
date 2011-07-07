/*
 *  vector_polynomial_gpu.h
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 18.03.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */

#ifndef VLI_VECTOR_POLYNOME_GPU_H
#define VLI_VECTOR_POLYNOME_GPU_H
#include "boost/swap.hpp"
#include "monome/polynome.h"
#include "monome/polynome_gpu.h"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "function_hooks/vli_vector_polynomial_gpu_function_hooks.hpp"

#include <ostream>


namespace vli
{

    // C I improved the inheritance, derived with one number and then realloc
    // C should we 100 % compatible with the stl::vector ?
    // C very bad news the size is static ................ 
    // C Bad luck due to my trick at the inialization I can not used the mother class ...
    // C Do you have any idea ?????
    // C Can we retemplate during run time compilation
    // C at the end I will think template is fortran programming .....
    // C But iw we consider the spare ou dense matrix on GPU with a fix size no pb ...... 
    // C I am preparign the worst case, mean dynamic

    template<class polynomial_gpu> 
    class vector_polynomial_gpu : public vli_gpu<typename polynomial_gpu::vli_value_type, 1>{ // 1 ....
    private:
        typedef typename polynomial_gpu::vli_value_type vli_value_type; // Just for convenience inside this class
        typedef typename std::size_t size_t;
        enum {max_order_poly = polynomial_gpu::max_order };
        enum {vli_size   = polynomial_gpu::size }; // C bad solution, I need the size of vli (#entry), because the polynomial has the good size
    public:
        // the usual proxy for have acces to element, we initialize with a polynomial, take care of the size !
        template<int OffSet>
        class proxy
        {
        public:
            proxy(vli_value_type* p, int i)
            :pdata_(p), pos(i), offset_(OffSet)
            {
            }
            
            proxy& operator=(polynomial_gpu const& p)
            {
                gpu::cu_check_error(cudaMemcpy((void*)(pdata_+pos*OffSet),(void*)(p.p()), max_order_poly*max_order_poly*vli_size*sizeof(typename polynomial_gpu::vli_value_type), cudaMemcpyDeviceToDevice ), __LINE__); 	           
                return *this;
            }
            
            friend std::ostream& operator << (std::ostream& os, proxy const & pr){
                pr.print(os);
                return os; 
            }
            
            void print(std::ostream & os) const{ // CONST ! 
                polynomial<vli_cpu<vli_value_type, vli_size>, max_order_poly > P;
                gpu::cu_check_error(cudaMemcpy((void*)(&P(0,0)),(void*)(pdata_+pos*OffSet),max_order_poly*max_order_poly*vli_size*sizeof(typename polynomial_gpu::vli_value_type), cudaMemcpyDeviceToHost), __LINE__);
                os << P;
            }
                                 
        private:
            int offset_; // it corresponds to the shift between each poly (arithmetic of pointer)
            vli_value_type* pdata_;
            int pos;
        };    
        
        vector_polynomial_gpu(size_t size)
        :size_(size),full_size_(size*max_order_poly*max_order_poly*vli_size)
        {
            this->realloc_lost(full_size_); // polynomial size : Order*Order*Vli::size
        }
        
        vector_polynomial_gpu()
        {
        //    this->realloc_lost(8); // crash if 0
        }

        vector_polynomial_gpu(vector_polynomial_gpu const& v) 
        {
            assert(v.size() == this->size()); 
            gpu::cu_check_error(cudaMemcpy((void*)this->p(), v.p(), full_size_*sizeof(typename polynomial_gpu::vli_value_type) , cudaMemcpyDeviceToDevice), __LINE__);
        }

        vector_polynomial_gpu& operator=(vector_polynomial_gpu v)
        {  
            swap(*this, v);
            return *this;
        }
        
        friend void swap(vector_polynomial_gpu & v1, vector_polynomial_gpu & v2)
        {
            using std::swap;
            swap(v1.data_, v2.data_);           
            swap(v1.size_,v2.size_);
            swap(v1.full_size_,v2.full_size_);
        }
        
        
        proxy<max_order_poly*max_order_poly*vli_size> operator[](size_t i)
        {
            return proxy<max_order_poly*max_order_poly*vli_size>(this->p(),i);
        }
    
        void resize(size_t num){ // bug if resize, TO DO !
            size_ = num;
            full_size_ = num*max_order_poly*max_order_poly*vli_size;
            this->realloc(full_size_); // if realloc smaller, data lost, done into vli_gpu
        } 
            
        size_t const & size() const {
             return size_;
        }
        
    private:
        size_t size_; // num of poly
        size_t full_size_; // number of VLI inside the vector
    };
    
    template <class BaseInt, int Size, int Order>
    vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> > inner_product(
                                                                                           vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> >  const& a, 
                                                                                           vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> >   const& b){
        assert(a.size() == b.size());
        vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> >  res;
        res.resize(a.size());
        inner_product_gpu(a,b,res);
        return res;
    }
    
    // C - tricky of tricky, is it the definition of french or bad programming ?
    template <class BaseInt, int Order, int SizeVli >
	std::ostream & operator<<(std::ostream & os, vector_polynomial_gpu< polynomial_gpu< vli_gpu<BaseInt, SizeVli>, Order > >   & v)
    {
        for(std::size_t i = 0; i < v.size(); i++)
            os << v[i] << std::endl;

        return os;
    }
}

#endif //VLI_VECTOR_POLYNOME_GPU_H