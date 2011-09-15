/*
 *  vector_polynomial_gpu.h
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 18.03.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */

#ifndef VLI_VECTOR_POLYNOMIAL_GPU_HPP
#define VLI_VECTOR_POLYNOMIAL_GPU_HPP
#include "vli/polynomial/polynomial_cpu.hpp"
#include "vli/polynomial/polynomial_gpu.hpp"
#include "vli/polynomial/monomial.hpp"
#include "vli/vli_cpu.hpp"
#include "vli/vli_gpu.hpp"
#include "vli/function_hooks/vli_vector_polynomial_gpu_function_hooks.hpp"
#include "vli/function_hooks/vli_number_gpu_function_hooks.hpp"
#include <boost/swap.hpp>
#include <ostream>





namespace vli
{

    template<class BaseInt, std::size_t Size>
    class vli_cpu;
    
    template<class vli_cpu, unsigned int Order>
	class polynomial_cpu;

	template <class vli_cpu>
    struct monomial;

    template<class polynomial>
    class vector_polynomial_cpu;

    template<class polynomial_gpu> 
    class vector_polynomial_gpu : public gpu_vector<typename polynomial_gpu::vli_value_type>{ 
    private:
        typedef typename polynomial_gpu::value_type::value_type vli_value_type; 
        enum {max_order_poly = polynomial_gpu::max_order };
        enum {vli_size   = polynomial_gpu::value_type::size }; 
        enum {element_offset = max_order_poly*max_order_poly*vli_size}; 
    public:
        typedef typename std::size_t size_type;

        // the usual proxy for have acces to element, we initialize with a polynomial, take care of the size !
        class proxy
        {
        public:
            proxy(vli_value_type* p, size_type i)
            :pdata_(p), pos_(i)
            {
            }
            
            operator polynomial_gpu() const // in case [] is right value
            {
                return polynomial_gpu(pdata_+pos_*element_offset);
            }
            
            proxy& operator=(polynomial_gpu const& p)
            {               
                gpu::cu_check_error(cudaMemcpy((void*)(pdata_+pos_*element_offset),(void*)p.p(),element_offset*sizeof(vli_value_type), cudaMemcpyDeviceToDevice ), __LINE__); 	           
                return *this;
            }
            
            proxy& operator=(proxy const& pr) // in case [] is right value
            {
                  gpu::cu_check_error(cudaMemcpy((void*)(pdata_+pos_*element_offset),(void*)(pr.pdata_+pr.pos_*element_offset),element_offset*sizeof(vli_value_type), cudaMemcpyDeviceToDevice ), __LINE__); 	           
                  return *this;
            }            
                       
            polynomial_gpu operator*(monomial<vli_gpu<vli_value_type,vli_size > > const& m) 
            {
                polynomial_gpu inter(pdata_+pos_*element_offset);
                inter *= m;  
                return inter;
            }
            
            polynomial_gpu operator*(int i)
            {
                vli_gpu<vli_value_type,vli_size >  inter(i);
                polynomial_gpu inter_poly(pdata_+pos_*element_offset);                
                inter_poly *= inter;
                return inter;
            }
            
            proxy& operator+=(polynomial_gpu const& p) 
            {
                //    plus_assign(this->p()+i*OffSet, m.coeff_.p());
                polynomial_gpu inter(pdata_+pos_*element_offset);
                inter += p; // I think by cheating I could merge this 2 lines
                *this = inter;
                return *this;
            }

            friend std::ostream& operator << (std::ostream& os, proxy const & pr){
                pr.print(os);
                return os; 
            }
            
            void print(std::ostream & os) const{ // CONST ! 
                polynomial_cpu<vli_cpu<vli_value_type, vli_size>, max_order_poly > P;
                gpu::cu_check_error(cudaMemcpy((void*)(&P(0,0)),(void*)(pdata_+pos_*element_offset),max_order_poly*max_order_poly*vli_size*sizeof(vli_value_type), cudaMemcpyDeviceToHost), __LINE__);
                os << P;
            }
            
            vli_gpu<vli_value_type, vli_size> BuildProxyToVli() const
            {
                vli_gpu<vli_value_type, vli_size> res;
                gpu::cu_check_error(cudaMemcpy((void*)(res.p()),(void*)(pdata_+pos_*element_offset), vli_size*sizeof(vli_value_type), cudaMemcpyDeviceToDevice ), __LINE__); 	           
                return res;
            }
            
                                 
        private:
            vli_value_type* pdata_;
            size_type pos_;
        };    
        
        vector_polynomial_gpu(size_type size = 1)
        : gpu_vector<typename polynomial_gpu::vli_value_type>(size*element_offset),size_(size)
        {   
        }
        
        vector_polynomial_gpu(vector_polynomial_gpu const& v)
        : gpu_vector<typename polynomial_gpu::vli_value_type>(v),size_(v.size_){	
           
        }
        
        /** CPU vector to GPU vector */
        explicit vector_polynomial_gpu(vector_polynomial_cpu< polynomial_cpu< vli_cpu <vli_value_type, vli_size>, max_order_poly > > const& vector){ 
            resize(vector.size()); // because default value is one !
            std::cout <<  vector.size()*max_order_poly*max_order_poly*vli_size*sizeof(vli_value_type) << std::endl;

            gpu::cu_check_error(cudaMemcpy( (void*)this->p(), (void*)&vector[0], vector.size()*max_order_poly*max_order_poly*vli_size*sizeof(vli_value_type), cudaMemcpyHostToDevice), __LINE__); 
        }

        vector_polynomial_gpu& operator=(vector_polynomial_gpu v)
        {
            swap(*this, v);
            return *this;
        }
                
        operator vector_polynomial_cpu<polynomial_cpu< vli_cpu<vli_value_type, vli_size>, max_order_poly> > () const
        {
            vector_polynomial_cpu< polynomial_cpu < vli_cpu<vli_value_type, vli_size>, max_order_poly> >  r(this->size());
            copy_vec_vli_to_cpu(r); //we are cheating 
            return r;
        }
        
        void copy_vec_vli_to_cpu( vector_polynomial_cpu< polynomial_cpu < vli_cpu<vli_value_type, vli_size>, max_order_poly> > & v) const
        {
            gpu::cu_check_error(cudaMemcpy( (void*)&v[0], (void*)this->p(),v.size()*element_offset*sizeof(vli_value_type),cudaMemcpyDeviceToHost ), __LINE__);					
        }
        
        void copy_vec_vli_to_gpu( vector_polynomial_cpu< polynomial_cpu < vli_cpu<vli_value_type, vli_size>, max_order_poly> > const & v) 
        {
            gpu::cu_check_error(cudaMemcpy( (void*)this->p(), (void*)&v[0],this->size()*element_offset*sizeof(vli_value_type),cudaMemcpyHostToDevice ), __LINE__);					
        }    
        
        proxy operator[](size_type i) 
        {
           return proxy(this->p(),i);
        }
        
        const proxy operator[](size_type i) const
        {
           return proxy(const_cast<vli_value_type*>(this->p()),i);
        }        
        
        inline size_type size() const{
            return size_;
        }
    
        void resize(size_type num){
            size_ = num;
            gpu_vector<typename polynomial_gpu::vli_value_type>::vec_resize(num*max_order_poly*max_order_poly*vli_size);
        }
        
        bool operator==(vector_polynomial_cpu<polynomial_cpu<vli_cpu<vli_value_type, vli_size>, max_order_poly> > const & v) const
        {
            return vector_polynomial_cpu<polynomial_cpu<vli_cpu<vli_value_type, vli_size>, max_order_poly> >(*this) == v;
        }

    private:
        size_type size_; // number of polynomial
        
    };
   
    template <class BaseInt, std::size_t Size, unsigned int Order>
    const polynomial_gpu<vli_gpu<BaseInt, Size>, Order>  
    inner_product( vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> >  const& a, 
                   vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> >  const& b){
        assert(a.size() == b.size());
        vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> > tmp;
        polynomial_gpu<vli_gpu<BaseInt, Size>, Order> res;
        inner_product_multiplication_gpu(a,b,tmp,res);
        return res;
    }
    
    template <class BaseInt, std::size_t Size, unsigned int Order >
	std::ostream & operator<<(std::ostream & os, vector_polynomial_gpu< polynomial_gpu< vli_gpu<BaseInt, Size>, Order > >   const& v)
    {
        for(std::size_t i = 0; i < v.size(); ++i)
            os << v[i] << std::endl;
        
        return os;
    }
}

#endif //VLI_VECTOR_POLYNOMIAL_GPU_HPP
