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
#include "polynomial/polynomial_cpu.hpp"
#include "polynomial/polynomial_gpu.hpp"
#include "polynomial/monomial.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "function_hooks/vli_vector_polynomial_gpu_function_hooks.hpp"
#include "function_hooks/vli_number_gpu_function_hooks.hpp"

#include <ostream>





namespace vli
{

    template<class BaseInt, int Size>
    class vli_cpu;
    
    template<class Vli, int Order>
	class polynomial_cpu;

	template <class Vli>
    struct monomial;

    template<class polynomial>
    class vector_polynomial_cpu;

    template<class polynomial_gpu> 
    class vector_polynomial_gpu : public gpu_vector<typename polynomial_gpu::vli_value_type>{ 
    private:
        typedef typename polynomial_gpu::vli_value_type vli_value_type; 
        typedef typename std::size_t size_t;
        enum {max_order_poly = polynomial_gpu::max_order };
        enum {vli_size   = polynomial_gpu::size }; 
        enum {OffSet = max_order_poly*max_order_poly*vli_size}; 
    public:
        typedef typename std::size_t size_type; // TO BE COMPLIANT WITH YOUR CODE ANDREAS
        // the usual proxy for have acces to element, we initialize with a polynomial, take care of the size !
        class proxy
        {
        public:
            proxy(vli_value_type* p, size_t i)
            :pdata_(p), pos(i)
            {
            }
            
            proxy& operator=(polynomial_gpu const& p)
            {
                gpu::cu_check_error(cudaMemcpy((void*)(pdata_+pos*OffSet),(void*)p.p(), max_order_poly*max_order_poly*vli_size*sizeof(typename polynomial_gpu::vli_value_type), cudaMemcpyDeviceToDevice ), __LINE__); 	           
                return *this;
            }
           
            template <class Vli> // Funny !
            proxy& operator+=(monomial<Vli> const& m)
            {
                assert(false);
            //    plus_assign(pdata_, m.coeff_.p());
                return *this;
            }
   
            friend std::ostream& operator << (std::ostream& os, proxy const & pr){
                pr.print(os);
                return os; 
            }
            
            void print(std::ostream & os) const{ // CONST ! 
                polynomial_cpu<vli_cpu<vli_value_type, vli_size>, max_order_poly > P;
                gpu::cu_check_error(cudaMemcpy((void*)(&P(0,0)),(void*)(pdata_+pos*OffSet),max_order_poly*max_order_poly*vli_size*sizeof(typename polynomial_gpu::vli_value_type), cudaMemcpyDeviceToHost), __LINE__);
                os << P;
            }
            
            vli_gpu<vli_value_type, vli_size> BuildProxyToVli() const
            {
                vli_gpu<vli_value_type, vli_size> res;
                gpu::cu_check_error(cudaMemcpy((void*)(res.p()),(void*)(pdata_+pos*OffSet), vli_size*sizeof(typename polynomial_gpu::vli_value_type), cudaMemcpyDeviceToDevice ), __LINE__); 	           
                return res;
            }
            
                                 
        private:
            vli_value_type* pdata_;
            size_t pos;
        };    
        
        vector_polynomial_gpu(size_t size = 1):gpu_vector<typename polynomial_gpu::vli_value_type>(size*OffSet),size_(size)
        {   
        }
        
        vector_polynomial_gpu(vector_polynomial_gpu const& v)
        :gpu_vector<typename polynomial_gpu::vli_value_type>(v),size_(v.size_){	
           
        }
        
        /** CPU vector to GPU vector */
        explicit vector_polynomial_gpu(vector_polynomial_cpu< polynomial_cpu< vli_cpu <vli_value_type, vli_size>, max_order_poly > >& vector){ 
             resize(vector.size()); // because default value is one !
             gpu::cu_check_error(cudaMemcpy( (void*)this->p(), (void*)&(vector[0](0,0)), vector.size()*max_order_poly*max_order_poly*vli_size*sizeof(vli_value_type), cudaMemcpyHostToDevice), __LINE__); 
        }

        vector_polynomial_gpu& operator=(vector_polynomial_gpu v)
        {  
            swap(*this, v);
            return *this;
        }
        
        operator vector_polynomial_cpu<polynomial_cpu< vli_cpu<vli_value_type, vli_size>, max_order_poly> > () const
        {
            vector_polynomial_cpu< polynomial_cpu < vli_cpu<vli_value_type, vli_size>, max_order_poly> >  r;
            copy_vec_vli_to_cpu(r);
            return r;
        }
        
        void copy_vec_vli_to_cpu( vector_polynomial_cpu< polynomial_cpu < vli_cpu<vli_value_type, vli_size>, max_order_poly> > & v) const
        {
            gpu::cu_check_error(cudaMemcpy( (void*)&v[0], (void*)this->p(),OffSet*sizeof(vli_value_type),cudaMemcpyDeviceToHost ), __LINE__);					
        }
        
        proxy operator[](size_t i)
        {
            return proxy(this->p(),i);
        }
    
        const proxy operator[](size_t i) const
        {
            return proxy(this->p(),i);
        }
        
        
        const size_t size() const{
            return size_;
        }
    
        void resize(std::size_t num){
            size_ = num;
            gpu_vector<typename polynomial_gpu::vli_value_type>::vec_resize(num*max_order_poly*max_order_poly*vli_size);
        }
        
        bool operator==(vector_polynomial_cpu<polynomial_cpu<vli_cpu<vli_value_type, vli_size>, max_order_poly> > const & v) const
        {
            return vector_polynomial_cpu<polynomial_cpu<vli_cpu<vli_value_type, vli_size>, max_order_poly> >(*this) == v;
        }
                                    
    private:
        std::size_t size_; // number of polynomial
        
    };
   
    template <class BaseInt, int Size, int Order>
    polynomial_gpu<vli_gpu<BaseInt, Size>, Order>  
    inner_product( vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> >  const& a, 
                   vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> >  const& b){
        assert(a.size() == b.size());
        polynomial_gpu<vli_gpu<BaseInt, Size>, Order> res;
        inner_product_multiplication_gpu(a,b,res);
        return res;
    }
    
    template <class BaseInt, int Order, int SizeVli >
	std::ostream & operator<<(std::ostream & os, vector_polynomial_gpu< polynomial_gpu< vli_gpu<BaseInt, SizeVli>, Order > >   const& v)
    {
        os << "---------- GPU ----------" << std::endl;
    
        for(std::size_t i = 0; i < v.size(); i++)
            os << v[i] << std::endl;

        return os;
    }
}

#endif //VLI_VECTOR_POLYNOME_GPU_H
