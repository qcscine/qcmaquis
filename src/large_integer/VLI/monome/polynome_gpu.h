/*
 *  monome_gpu.h
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 18.03.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */

#ifndef VLI_POLYNOME_GPU_H
#define VLI_POLYNOME_GPU_H
#include "boost/swap.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "function_hooks/vli_polynomial_gpu_function_hooks.hpp"
#include "monome/monome.h"

#include <ostream>
#include <cmath>

namespace vli
{	    
    template<class BaseInt, int Size>
    class vli_cpu;
    
    
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
    
    /**
     * Multiplication of  polynomial * monomial
     */
    template<class Vli, int Order>
    polynomial_gpu<Vli, Order> operator * (polynomial_gpu<Vli, Order> const& p1, monomial<Vli> const &m2)
    {
        polynomial_gpu<Vli, Order> result;
        poly_mono_multiply(result, p1, m2);
        return result;
    }
    
    
	template<class Vli, int Order>
	class polynomial_gpu : public vli_gpu<typename Vli::value_type, Order*Order*Vli::size> {
        public:
 
        typedef typename Vli::value_type vli_value_type; // Just for convenience inside this class
        typedef typename Vli::size_type size_type; // Type of the exponents (has to be the same type as Vli::size_type)
        enum { max_order = Order };
        enum { size = Vli::size }; // C to remove if you find a solution see vector_poly ...

        class proxy
        {
        public:
            proxy(polynomial_gpu& poly, size_type j, size_type h)
            :data_(poly.p()),j_(j),h_(h){}
            
            proxy& operator= (Vli const& vli ){
                gpu::cu_check_error(cudaMemcpy((void*)(data_+(j_*max_order*Vli::size+h_*Vli::size)), (void*)vli.p(), Vli::size*sizeof(vli_value_type), cudaMemcpyDeviceToDevice), __LINE__);
                return *this;
            }
       
           //TODO do we really want this? Or just using an explicit conversion? 
           // C - you have to explain me 
            proxy& operator= (vli_cpu<vli_value_type, Vli::size> const& vli ){
                gpu::cu_check_error(cudaMemcpy((void*)(data_+(j_*max_order*Vli::size+h_*Vli::size)), (void*)&vli[0], Vli::size*sizeof(vli_value_type), cudaMemcpyHostToDevice), __LINE__);
                return *this;
            }

            typename vli_gpu<typename Vli::value_type,Order*Order*Vli::size>::proxy operator[](size_type i)
            {
                return typename vli_gpu<typename Vli::value_type,Order*Order*Vli::size>::proxy(data_+(j_*max_order*Vli::size+h_*Vli::size),i);
            }

            friend std::ostream& operator << (std::ostream& os, proxy const& pr){
                pr.print(os);
                return os;
            }
            
            void print(std::ostream& os) const{
                vli_cpu<vli_value_type,Vli::size> vli;
                gpu::cu_check_error(cudaMemcpy((void*)&vli[0],(void*)(data_+(j_*max_order*Vli::size+h_*Vli::size)),Vli::size*sizeof(vli_value_type), cudaMemcpyDeviceToHost), __LINE__);
                os << vli;
            }
            
        private:
            vli_value_type* data_;   
            size_type j_; // for operator (j,h) 
            size_type h_; // for operator (j,h)
        };
        /** Serious on the specialization */
        friend void poly_multiply <>(polynomial_gpu& result , polynomial_gpu const& p1, polynomial_gpu const& p2);
        friend void poly_mono_multiply <>(polynomial_gpu& result , polynomial_gpu const& p1, monomial<Vli>  const& m2);
       
        
        polynomial_gpu(){
        }
        
        /** GPU poly to GPU poly */
        polynomial_gpu(polynomial<vli_cpu<vli_value_type, Vli::size>, Order>& poly){ 
            gpu::cu_check_error(cudaMemcpy( (void*)this->p(), (void*)&poly(0,0), Order*Order*Vli::size*sizeof(vli_value_type), cudaMemcpyHostToDevice), __LINE__); 
        }
        
        operator polynomial<vli_cpu<vli_value_type, Vli::size>, Order>() const
        {
            polynomial<vli_cpu<vli_value_type, Vli::size>, Order> r;
            copy_poly_vli_to_cpu(r);
            return r;
        }
        
        void copy_poly_vli_to_cpu(polynomial<vli_cpu<vli_value_type, Vli::size>, Order> & p) const
        {
            gpu::cu_check_error(cudaMemcpy( (void*)&p(0,0)[0], (void*)this->p(), Order*Order*Vli::size*sizeof(vli_value_type),cudaMemcpyDeviceToHost ), __LINE__);					
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
        
        /**  TO DO shame on me, think to integrate *= += into the kernel **/
        polynomial_gpu operator *= (Vli const& c)
        {
            polynomial_gpu<Vli, Order> result;
            poly_mono_multiply(result, (*this), c);
            *this = result;
            return *this;
        }
        
        polynomial_gpu& operator *= (monomial<Vli> const& c)
        {
            (*this) *= c.coeff;
            return *this;
        }
        
        /** GPU/GPU**/
        bool operator==(polynomial_gpu const & p) const 
        {
            //TODO
            //BUGY , presently the kernel destroy the data, I do not know why !
            //assert(false);
            int test(0); // shoud be a bool
            detail::equal_gpu((*this).p(), p.p(), max_order*max_order*Vli::size, &test);    
            return (test == 0) ? true : false;
        }

         /** GPU/CPU, order cares !**/
        bool operator==(polynomial<vli_cpu<vli_value_type, Vli::size>, Order>  & p) const
        {
            return polynomial<vli_cpu<vli_value_type, Vli::size>, Order > (*this) == p;
         //   return (*this) == polynomial_gpu(p);
        }
     
        proxy operator ()(size_type j_exp, size_type h_exp) 
        {
            assert(j_exp < max_order);
            assert(h_exp < max_order);
            return proxy(*this, j_exp, h_exp);
        }
        
        /** 
         Due to the derivation a large part of the operators come from vli_gpu
         */
    }; //end class
    
    /**
     * Multiplication of a polynomial with a factor
     */
    template<class Vli, int Order>
    polynomial_gpu<Vli, Order> operator * (polynomial_gpu<Vli, Order> p, Vli const& c)
    {
        p *= c;
        return p;
    }
    
    template<class Vli, int Order>
    polynomial_gpu<Vli, Order> operator * (Vli const& c, polynomial_gpu<Vli, Order> const& p)
    {
        return p * c;
    }
    
    template <class BaseInt, int Order, int SizeVli >
	std::ostream & operator<<(std::ostream & os, polynomial_gpu< vli_gpu<BaseInt, SizeVli>, Order > & p)
    {
        os << "------------------------" << std::endl;
        for(std::size_t i = 0; i < Order ; i++){
            for(std::size_t j = 0; j < Order ; j++){
                os << "gpu --- Coeff (j,h) = " << i <<" "<<j<< std::endl;
                os << p(i,j) << std::endl;
            }
        }
        os <<  "------------------------" << std::endl;
                
        return os;
    }
    

}

#endif //VLI_MONOME_GPU_H
