#ifndef INNER_PRODUCT_GPU_BOOSTER_HPP
#define INNER_PRODUCT_GPU_BOOSTER_HPP

#include "vli/utils/gpu_manager.h"
#include "utils/timings.h"

#include "vli/detail/kernels_gpu.h"
#include <iostream>

namespace vli
{

template <class BaseInt, std::size_t Size>
class vli_cpu;

template <class Vli, unsigned int Order>
class polynomial;

template <class Polynomial>
class vector_polynomial;

namespace detail
{
    template <typename T>
    class gpu_memblock
    {
      public:
        typedef T value_type;
        gpu_memblock(std::size_t size) : size_(size), data_(NULL) {
            gpu::cu_check_error(cudaMalloc((void**)&(this->data_), size*sizeof(T)), __LINE__);
        }

        ~gpu_memblock() {
            cudaFree(this->data_);
        }

        std::size_t size() const {
            return size_;
        }

        T* p() {
            return data_;
        }

        T const* p() const {
            return data_;
        }

      private:
        std::size_t size_;
        T* data_;
    };
    
    template <typename Vli>
    class asm_gpu_add{
      private:
        typedef typename Vli::value_type  base_int_type;
        enum {Size =Vli::size};
        enum {Order =11};
      public:
        explicit asm_gpu_add(vli_cpu<base_int_type, Size> &a, vli_cpu<base_int_type, Size> const &b):v1_(Size),v2_(Size){
            gpu::cu_check_error(cudaMemcpyAsync((void*)v1_.p(),(void*)&a[0],Size*sizeof(base_int_type),cudaMemcpyHostToDevice),__LINE__);
            gpu::cu_check_error(cudaMemcpyAsync((void*)v2_.p(),(void*)&b[0],Size*sizeof(base_int_type),cudaMemcpyHostToDevice),__LINE__);
            addition_gpu(vli_size_tag<Vli::size,Order>(), v1_.p(), v2_.p());
            gpu::cu_check_error(cudaGetLastError(),__LINE__);
        } 
        operator vli_cpu<typename Vli::value_type, Vli::size >() const {
            vli_cpu<typename Vli::value_type, Vli::size> res;
            gpu::cu_check_error(cudaMemcpy((void*)&res[0],(void*)v1_.p(),Size*sizeof(base_int_type),cudaMemcpyDeviceToHost),__LINE__);
            return res;
        }
      private:
        gpu_memblock<base_int_type> v1_;
        gpu_memblock<base_int_type> v2_;
    };

    template <typename Vli>
    class asm_gpu_mul{
      private:
        typedef typename Vli::value_type base_int_type;
        enum {Size = Vli::size};
        enum {Order =11};
      public:
        explicit asm_gpu_mul(vli_cpu<base_int_type, 2*Size> &a, vli_cpu<base_int_type, Size> const &b, vli_cpu<base_int_type, Size> const &c):v1_(2*Size),v2_(Size),v3_(Size){
            gpu::cu_check_error(cudaMemcpyAsync((void*)v1_.p(),(void*)&a[0],2*Size*sizeof(base_int_type),cudaMemcpyHostToDevice),__LINE__);
            gpu::cu_check_error(cudaMemcpyAsync((void*)v2_.p(),(void*)&b[0],Size*sizeof(base_int_type),cudaMemcpyHostToDevice),__LINE__);
            gpu::cu_check_error(cudaMemcpyAsync((void*)v3_.p(),(void*)&c[0],Size*sizeof(base_int_type),cudaMemcpyHostToDevice),__LINE__);
            multiplication_gpu(vli_size_tag<Size,Order>(), v1_.p(), v2_.p(), v3_.p());
            gpu::cu_check_error(cudaGetLastError(),__LINE__);
        } 
        operator vli_cpu<base_int_type, 2*Size >() const {
            vli_cpu<base_int_type, 2*Size> res;
            gpu::cu_check_error(cudaMemcpy((void*)&res[0],(void*)v1_.p(),2*Size*sizeof(base_int_type),cudaMemcpyDeviceToHost),__LINE__);
            return res;
        }
      private:
        gpu_memblock<base_int_type> v1_;
        gpu_memblock<base_int_type> v2_;
        gpu_memblock<base_int_type> v3_;
    };


    // TODO this template parameters should be more restrictive
    template <typename Vli, unsigned int Order>
    class inner_product_gpu_booster
    {
      private:
        typedef typename Vli::value_type           base_int_type;
        enum {factor_element_size = Order * Order * Vli::size };
        enum {product_element_size = 2*Order * 2*Order * 2 * Vli::size }; // VLi are twice larger
      public:
        typedef std::size_t size_type;
   
        explicit inner_product_gpu_booster(
                  vector_polynomial<polynomial<Vli,Order> > const& v1
                , vector_polynomial<polynomial<Vli,Order> > const& v2
                , size_type partsize
                )
        : v1_(partsize*factor_element_size), v2_(partsize*factor_element_size), tmp_(partsize*product_element_size), partsize_(partsize)
        {
            assert(partsize <= v1.size());
            assert(partsize <= v2.size());

            gpu::cu_check_error(cudaMemcpyAsync((void*)v1_.p(),(void*)&v1[0],partsize*factor_element_size*sizeof(base_int_type),cudaMemcpyHostToDevice),__LINE__);
            gpu::cu_check_error(cudaMemcpyAsync((void*)v2_.p(),(void*)&v2[0],partsize*factor_element_size*sizeof(base_int_type),cudaMemcpyHostToDevice),__LINE__);
            gpu::cu_check_error(cudaMemset((void*)tmp_.p(),0,partsize*product_element_size*sizeof(base_int_type)),__LINE__);

            // run/queue the inner product computation
            inner_product_vector(vli_size_tag<Vli::size,Order>(), partsize, v1_.p(), v2_.p(), tmp_.p());
            gpu::cu_check_error(cudaGetLastError(),__LINE__);
        }
        
        operator polynomial<vli_cpu<typename Vli::value_type,  2*Vli::size >,2*Order>() const {
            polynomial< vli_cpu<typename Vli::value_type,  2*Vli::size >,2*Order> poly; 
            gpu::cu_check_error(cudaMemcpy((void*)&poly(0,0),(void*)tmp_.p(),product_element_size*sizeof(base_int_type),cudaMemcpyDeviceToHost),__LINE__);
            return poly;
        }

      private:
        gpu_memblock<base_int_type> v1_;
        gpu_memblock<base_int_type> v2_;
        gpu_memblock<base_int_type> tmp_;
        std::size_t partsize_;
    };
}

namespace detail
{
#ifdef _OPENMP
template <class BaseInt, std::size_t Size, unsigned int Order>
polynomial<vli_cpu<BaseInt, 2*Size>, 2*Order> 
inner_product_openmp_gpu( vector_polynomial<polynomial<vli_cpu<BaseInt, Size>, Order> >  const& v1, 
               vector_polynomial<polynomial<vli_cpu<BaseInt, Size>, Order> >  const& v2){

    std::cout<<"inner_product: OpenMP +CUDA"<<std::endl;
    assert(v1.size() == v2.size());
    std::size_t size_v = v1.size();
    polynomial<vli_cpu<BaseInt, 2*Size>, 2*Order>  res[omp_get_max_threads()];
    std::size_t split = (std::size_t)(v1.size()*1);
    detail::inner_product_gpu_booster<vli_cpu<BaseInt,Size>,Order> gpu_product(v1,v2,split);

    #pragma omp parallel for
    for(std::size_t i=split ; i < size_v ; ++i)
        res[omp_get_thread_num()] += v1[i]*v2[i];

    for(std::size_t i=1; i < omp_get_max_threads(); ++i)
        res[0]+=res[i];

    res[0] += polynomial<vli_cpu<BaseInt, 2*Size>, 2*Order >(gpu_product); // this thing synchronize
    return res[0];
}
#endif

template <class BaseInt, std::size_t Size, unsigned int Order>
polynomial<vli_cpu<BaseInt, 2*Size>, 2*Order> 
inner_product_gpu( vector_polynomial<polynomial<vli_cpu<BaseInt, Size>, Order> >  const& v1, 
               vector_polynomial<polynomial<vli_cpu<BaseInt, Size>, Order> >  const& v2){
    std::cout<<"inner_product: single thread + CUDA"<<std::endl;
    assert(v1.size() == v2.size());
    std::size_t size_v = v1.size();
    
    polynomial<vli_cpu<BaseInt,2*Size>, 2*Order> res;
    std::size_t split = static_cast<std::size_t>(1*v1.size());
     
    detail::inner_product_gpu_booster<vli_cpu<BaseInt,Size>,Order> gpu_product(v1,v2,split);

    for(std::size_t i=split ; i < size_v ; ++i){
        res += v1[i]*v2[i];
    }
   
    res += polynomial<vli_cpu<BaseInt, 2*Size>, 2*Order >(gpu_product);
    
    return res;
}

template<class BaseInt, std::size_t Size>
vli_cpu<BaseInt, Size> addition_gpu(vli_cpu<BaseInt, Size> & a, vli_cpu<BaseInt, Size> const&b){
    asm_gpu_add<vli_cpu<BaseInt,Size> > addition(a,b); 
    return vli_cpu<BaseInt, Size> (addition);
} 

template<class BaseInt, std::size_t Size>
vli_cpu<BaseInt, 2*Size> multiplication_gpu(vli_cpu<BaseInt,2*Size> & a, vli_cpu<BaseInt, Size> const&b, vli_cpu<BaseInt, Size> const&c){
    asm_gpu_mul<vli_cpu<BaseInt,Size> > multiplication(a,b,c);
    return vli_cpu<BaseInt, 2*Size> (multiplication);
} 

} // end namespace detail
} // end namespace vli

#endif //INNER_PRODUCT_GPU_BOOSTER_HPP
