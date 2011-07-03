/*
 *  monome.h
 *  vli
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 18.03.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */

#ifndef VLI_MONOME_GPU_H
#define VLI_MONOME_GPU_H
#include "boost/swap.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "detail/vli_polynomial_gpu_function_hooks.hpp"

#include <ostream>
#include <cmath>

namespace vli
{	    
    template<class BaseInt, int Size>
    class vli_cpu;
    
//    template<class BaseInt, int Size>
//    class vli_gpu;
    
    template<class Vli, int Order>
	class polynomial;
    
    /**
     * Multiplication of two polynomials
     */
   	template<class Vli, int Order>
    polynomial_gpu<Vli, Order> operator * (polynomial_gpu<Vli, Order> const& p1, polynomial_gpu<Vli, Order> const& p2)
    {
        polynomial_gpu<Vli, Order> result;
        poly_multiply(result, p1, p2);
        return result;
    }
    
	template<class Vli, int Order>
	class polynomial_gpu : public vli_gpu<typename Vli::value_type, Order*Order*Vli::vli_size> {
        public:
        typedef typename Vli::value_type value_type;
        typedef typename Vli::size_type size_type;	
        enum { vli_size  = Vli::vli_size, max_order = Order };
        
        class proxy
        {
        public:
            proxy(polynomial_gpu& poly, size_type j, size_type h)
            :data_(poly.p()),j_(j),h_(h),i_(-1){}
            
            proxy& operator= (vli_gpu<value_type, vli_size> & vli ){
                gpu::cu_check_error(cudaMemcpy((void*)(data_+(j_*max_order*vli_size+h_*vli_size)), (void*)vli.p(), vli_size*sizeof(value_type), cudaMemcpyDeviceToDevice), __LINE__);
                return *this;
            }
        
            proxy& operator= (vli_cpu<value_type, vli_size> & vli ){
                gpu::cu_check_error(cudaMemcpy((void*)(data_+(j_*max_order*vli_size+h_*vli_size)), (void*)&vli[0], vli_size*sizeof(value_type), cudaMemcpyHostToDevice), __LINE__);
                return *this;
            }
            
            proxy& operator= (value_type i){
                if(i_ != -1){
                    value_type num(i);
                gpu::cu_check_error(cudaMemcpy((void*)(data_+(j_*max_order*vli_size+h_*vli_size+i_)), (void*)&num, sizeof(value_type), cudaMemcpyHostToDevice), __LINE__);
                }else{
                    assert(false); //dummy case, e.g. pa(0,0) = 255; or set everything to the corresponding value
                }
                return *this;
            }
            
            proxy& operator[] (size_type i ){
                i_=i;
                return *this;
            }
                           
            friend std::ostream& operator << (std::ostream& os, proxy const& pr){
                pr.print(os);
                return os;
            }
            
            void print(std::ostream& os) const{
                vli_cpu<value_type,vli_size> vli;
                gpu::cu_check_error(cudaMemcpy((void*)&vli[0],(void*)(data_+(j_*max_order*vli_size+h_*vli_size)),vli_size*sizeof(value_type), cudaMemcpyDeviceToHost), __LINE__);
                os << vli;
            }
            
        private:
            value_type* data_;   
            size_type i_; // for operator [i]
            size_type j_; // for operator (j,h) 
            size_type h_; // for operator (j,h)
        };
        
        friend polynomial_gpu operator * <> (polynomial_gpu const& p, polynomial_gpu const& m);
        friend void poly_multiply <>(polynomial_gpu& result , polynomial_gpu const& p1, polynomial_gpu const& p2);
       
        polynomial_gpu(){}
        
        /** GPU poly to GPU poly */
        polynomial_gpu(polynomial<vli_cpu<value_type, vli_size>, Order>& poly){ 
            gpu::cu_check_error(cudaMemcpy( (void*)this->p(), (void*)&poly(0,0), Order*Order*vli_size*sizeof(value_type), cudaMemcpyHostToDevice), __LINE__); 
        }
        
        operator polynomial<vli_cpu<value_type, vli_size>, Order>() const
        {
            polynomial<vli_cpu<value_type, vli_size>, Order> r;
            copy_poly_vli_to_cpu(r);
            return r;
        }
        
        void copy_poly_vli_to_cpu(polynomial<vli_cpu<value_type, vli_size>, Order> & p) const
        {
            gpu::cu_check_error(cudaMemcpy( (void*)&p(0,0)[0], (void*)this->p(), Order*Order*vli_size*sizeof(value_type),cudaMemcpyDeviceToHost ), __LINE__);					
        }
        
        inline Vli operator ()(unsigned int j_exp, unsigned int h_exp) const
        {
            assert(j_exp < max_order);
            assert(h_exp < max_order);
            return (*this)[j_exp*max_order+h_exp];
        }
       
        
        /**
        * Plus assign with a polynomial
        */
        polynomial_gpu& operator += (polynomial_gpu const& p)
        {
            poly_addition(*(this), p);
            return *this;
        }
        
        
        /** GPU/GPU**/
        bool operator==(polynomial_gpu const & p) const 
        {
            //BUGY , presently the kernel destroy the data, I do not know why !
            assert(false);
            int test(0); // shoud be a bool
            detail::equal_gpu((*this).p(), p.p(), max_order*max_order*vli_size, &test);    
            return (test == 0) ? true : false;
        }

         /** GPU/CPU, order cares !**/
        bool operator==(polynomial<vli_cpu<value_type, vli_size>, Order>  & p) const
        {
            return polynomial<vli_cpu<value_type, vli_size>, Order > (*this) == p;
         //   return (*this) == polynomial_gpu(p);
        }
     
        proxy operator ()(size_type j_exp, size_type h_exp) 
        {
            assert(j_exp < max_order);
            assert(h_exp < max_order);
            return proxy(*this, j_exp, h_exp);
        }
        
        /** 
         Due to the derivation a large part of the operator come from vli_gpu
         */

    };
}

#endif //VLI_MONOME_GPU_H
