/*
 *  monome_gpu.h
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 18.03.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */

#ifndef VLI_POLYNOMIAL_GPU_HPP
#define VLI_POLYNOMIAL_GPU_HPP
#include "vli/vli_gpu.hpp"
#include "vli/function_hooks/vli_polynomial_gpu_function_hooks.hpp"
#include "vli/polynomial/monomial.hpp"

#include <boost/swap.hpp>
#include <ostream>
#include <cmath>

namespace vli
{	    
    template<class BaseInt, std::size_t Size>
    class vli_cpu;
    
    
    template<class Vli, unsigned int Order>
	class polynomial_cpu;
    
    
    
    /**
     * Multiplication of two polynomials
     */
   	template<class Vli, unsigned int Order>
    const polynomial_gpu<Vli, Order> operator * (polynomial_gpu<Vli, Order> const& p1, polynomial_gpu<Vli, Order> const& p2)
    {
        polynomial_gpu<Vli, Order> result;
        poly_multiply(result, p1, p2);
        return result;
    }
    
    /**
     * Multiplication of  polynomial_cpu * monomial
     */
    template<class Vli, unsigned int Order>
    const polynomial_gpu<Vli, Order> operator * (polynomial_gpu<Vli, Order> const& p1, monomial<Vli> const &m2)
    {
        polynomial_gpu<Vli, Order> result;
        poly_mono_multiply(result, p1, m2);
        return result;
    }
    
    
	template<class Vli, unsigned int Order>
	class polynomial_gpu :  public  gpu_array<typename Vli::value_type,  Order*Order*Vli::size> 
    {
        public:
 
        typedef typename Vli::value_type vli_value_type; // Just for convenience inside this class
        typedef unsigned int exponent_type; // Type of the exponents (has to be the same type as Vli::size_type)
        typedef Vli value_type;
        enum { max_order = Order };
//        enum { size = Vli::size }; // C to remove if you find a solution see vector_poly ...

        class proxy
        {
        public:
            proxy(polynomial_gpu& poly, exponent_type j, exponent_type h)
            :data_(poly.p()),j_(j),h_(h){}
            
            proxy& operator= (Vli const& vli ){
                gpu::cu_check_error(cudaMemcpy((void*)(data_+(j_*max_order*Vli::size+h_*Vli::size)), (void*)vli.p(), Vli::size*sizeof(vli_value_type), cudaMemcpyDeviceToDevice), __LINE__);
                return *this;
            }
            
            bool operator>(int i){
                vli_cpu<vli_value_type,Vli::size> vli;
                gpu::cu_check_error(cudaMemcpy((void*)&vli[0],(void*)(data_+(j_*max_order*Vli::size+h_*Vli::size)),Vli::size*sizeof(vli_value_type), cudaMemcpyDeviceToHost), __LINE__);
                return (vli > i);
            }

            bool operator<(int i){
                vli_cpu<vli_value_type,Vli::size> vli;
                gpu::cu_check_error(cudaMemcpy((void*)&vli[0],(void*)(data_+(j_*max_order*Vli::size+h_*Vli::size)),Vli::size*sizeof(vli_value_type), cudaMemcpyDeviceToHost), __LINE__);
                return (vli < i);
            }
       
            typename vli_gpu<typename Vli::value_type,Order*Order*Vli::size>::proxy operator[](typename Vli::size_type i)
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
            exponent_type j_; // for operator (j,h)
            exponent_type h_; // for operator (j,h)
        };
        /** Serious on the specialization */
      //friend void poly_multiply <>(polynomial_gpu& result , polynomial_gpu const& p1, polynomial_gpu const& p2);
      //   friend void poly_mono_multiply <>(polynomial_gpu& result , polynomial_gpu const& p1, monomial<Vli>  const& m2);
       
        
        polynomial_gpu(){
        }
        
        explicit polynomial_gpu(vli_value_type* p){
            gpu::cu_check_error(cudaMemcpy((void*)(this->p()), (void*)p, Order*Order*Vli::size*sizeof(vli_value_type), cudaMemcpyDeviceToDevice ), __LINE__);
        }
        
        /** CPU poly to GPU poly */
        explicit polynomial_gpu(polynomial_cpu<vli_cpu<vli_value_type, Vli::size>, Order> const& poly){ 
            gpu::cu_check_error(cudaMemcpy( (void*)this->p(), (void*)&(poly(0,0)[0]), Order*Order*Vli::size*sizeof(vli_value_type), cudaMemcpyHostToDevice), __LINE__);             
        }
        
        operator polynomial_cpu<vli_cpu<vli_value_type, Vli::size>, Order>() const
        {
            polynomial_cpu<vli_cpu<vli_value_type, Vli::size>, Order> r;
            copy_poly_vli_to_cpu(r);
            return r;
        }
        
        void copy_poly_vli_to_cpu(polynomial_cpu<vli_cpu<vli_value_type, Vli::size>, Order> & p) const
        {
            gpu::cu_check_error(cudaMemcpy( (void*)&p(0,0)[0], (void*)this->p(), Order*Order*Vli::size*sizeof(vli_value_type),cudaMemcpyDeviceToHost ), __LINE__);					
        }
               
        /**
        * Plus assign with a polynomial_cpu
        */
        polynomial_gpu& operator += (polynomial_gpu const& p)
        {
            plus_assign_poly(*(this), p);
            return *this;
        }        
        
        template <typename T>
        polynomial_gpu& operator += (T const& t)
        {
            // Possible because we sum only the first coefficient
            using vli::poly_addition_int;
            poly_addition_int(*this,t);
            return *this;
        }
        
        
        polynomial_gpu& operator += (monomial<Vli> const& m)
        {
            assert(false);
            return *this;
        }
        
        
        /**
         * Minus assign with a polynomial_cpu
         */
        polynomial_gpu& operator -= (polynomial_gpu const& p)
        {
            poly_substraction(*(this), p);
            return *this;
        }
        
        
        /**  TO DO shame on me, think to integrate *= += into the kernel **/
        polynomial_gpu& operator *= (Vli const& c)
        {
            polynomial_gpu<Vli, Order> result;
            monomial<Vli> m(c,0,0);
            poly_mono_multiply(result, (*this), m);
            *this = result;
            return *this;
        }
        
        polynomial_gpu& operator *= (monomial<Vli> const& m)
        {            // TO DO TO CHANGE
            polynomial_gpu<Vli, Order> result;
            poly_mono_multiply(result, *this, m);
            *this = result;
            return *this;
        }
        
         /** GPU/CPU, order cares !**/
        bool operator==(polynomial_cpu<vli_cpu<vli_value_type, Vli::size>, Order>  & p) const
        {
            return polynomial_cpu<vli_cpu<vli_value_type, Vli::size>, Order > (*this) == p;
        }
        
     
        proxy operator ()(exponent_type j_exp, exponent_type h_exp) 
        {
            assert(j_exp < max_order);
            assert(h_exp < max_order);
            return proxy(*this, j_exp, h_exp);
        }
        
        
        const Vli operator ()(exponent_type j_exp, exponent_type h_exp) const
        {
            assert(j_exp < max_order);
            assert(h_exp < max_order);
            return Vli(this->data_+(j_exp*max_order+h_exp)*Vli::size);
        }
        
        /** 
         Due to the derivation a large part of the operators come from vli_gpu
         */
    }; //end class
    
    /**
     * Multiplication of a polynomial_cpu with a factor
     */
    template<class Vli, unsigned int Order>
    polynomial_gpu<Vli, Order> operator * (polynomial_gpu<Vli, Order> p, Vli const& c)
    {
        p *= c;
        return p;
    }
    
    template<class Vli, unsigned int Order>
    polynomial_gpu<Vli, Order> operator * (Vli const& c, polynomial_gpu<Vli, Order> const& p)
    {
        return p * c;
    }
    
    template<class Vli, unsigned int Order>
    polynomial_gpu<Vli, Order> operator * (int c, polynomial_gpu<Vli, Order> const& p)
    {
        return Vli(c) * p;
    }
    
    
    template <class BaseInt, std::size_t Size, unsigned int Order>
	std::ostream & operator<<(std::ostream & os, polynomial_gpu< vli_gpu<BaseInt, Size>, Order > const& p)
    {
        typedef typename polynomial_gpu<vli_gpu<BaseInt,Size>,Order>::exponent_type exponent_type;
        for(exponent_type i = 0; i < Order ; ++i)
            for(exponent_type j = 0; j < Order ; ++j)
                os << p(i,j) << std::endl;
        return os;
    }
    

}

#endif //VLI_POLYNOMIAL_GPU_HPP
