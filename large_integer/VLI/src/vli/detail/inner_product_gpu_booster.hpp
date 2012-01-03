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
class polynomial_cpu;

template <class Polynomial>
class vector_polynomial_cpu;

namespace detail
{
    template <typename T>
    class gpu_memblock
    {
      public:
        typedef T value_type;
        gpu_memblock(std::size_t size)
            :size_(size), data_(NULL)
        {
            gpu::cu_check_error(cudaMalloc((void**)&(this->data_), size*sizeof(T)), __LINE__);
        }

        ~gpu_memblock()
        {
            cudaFree(this->data_);
        }

        std::size_t size() const
        {
            return size_;
        }

        T* p()
        {
            return data_;
        }

        T const* p() const
        {
            return data_;
        }

      private:
        std::size_t size_;
        T* data_;

    };
    
    // TODO this template parameters should be more restrictive
    template <typename Vli, unsigned int Order>
    class inner_product_gpu_booster
    {
      private:
//        typedef typename VectorPolynomial::value_type   polynomial_type;
//        typedef typename polynomial_type::value_type    vli_type;
        typedef typename Vli::value_type           base_int_type;
        enum {factor_element_size = Order * Order * Vli::size };
        enum {product_element_size = 2*Order * 2*Order * Vli::size };
        enum {num_threads = 64}; // TODO change
      public:
        typedef std::size_t size_type;
   
        explicit inner_product_gpu_booster(
                  vector_polynomial_cpu<polynomial_cpu<Vli,Order> > const& v1
                , vector_polynomial_cpu<polynomial_cpu<Vli,Order> > const& v2
                , size_type partsize
                )
        : v1_(partsize*factor_element_size), v2_(partsize*factor_element_size), tmp_(partsize*product_element_size), partsize_(partsize)
        {
            assert(partsize <= v1.size());
            assert(partsize <= v2.size());
            gpu::cu_check_error(cudaMemcpyAsync((void*)v1_.p(),(void*)&v1[0],partsize*factor_element_size*sizeof(base_int_type),cudaMemcpyHostToDevice), __LINE__);
            gpu::cu_check_error(cudaMemcpyAsync((void*)v2_.p(),(void*)&v2[0],partsize*factor_element_size*sizeof(base_int_type),cudaMemcpyHostToDevice),__LINE__);
            gpu::cu_check_error(cudaMemset((void*)tmp_.p(),0,partsize*product_element_size*sizeof(base_int_type)),__LINE__);

            // run/queue the inner product computation
            inner_product_vector(vli_size_tag<Vli::size>(), Order, partsize, v1_.p(), v2_.p(), tmp_.p(), num_threads);
            gpu::cu_check_error(cudaGetLastError(),__LINE__);

        }
        
        operator polynomial_cpu<Vli,2*Order>() const
        {
            std::vector<polynomial_cpu<Vli,2*Order> > restmp(partsize_); 
            polynomial_cpu<Vli,2*Order> poly; 

            gpu::cu_check_error(cudaMemcpy((void*)&restmp[0](0,0),(void*)tmp_.p(),partsize_*product_element_size*sizeof(base_int_type),cudaMemcpyDeviceToHost),__LINE__);
            // C - a reduction of 16384 poly  order 21*21 costs 0.1 s .... so CPU
            for(std::size_t i=0; i <partsize_;++i){
               poly += restmp[i];
            }
 
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
polynomial_cpu<vli_cpu<BaseInt, Size>, 2*Order> 
inner_product_openmp_gpu( vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v1, 
               vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v2){

    std::cout<<"inner_product: OpenMP +CUDA"<<std::endl;
    assert(v1.size() == v2.size());
    std::size_t size_v = v1.size();
    polynomial_cpu<vli_cpu<BaseInt, Size>, 2*Order>  res[omp_get_max_threads()];
    std::size_t split = static_cast<std::size_t>(v1.size()*VLI_SPLIT_PARAM);
    detail::inner_product_gpu_booster<vli_cpu<BaseInt,Size>,Order> gpu_product(v1,v2,split);

    #pragma omp parallel for
    for(std::size_t i=split ; i < size_v ; ++i){
        res[omp_get_thread_num()] += v1[i]*v2[i];
    }

    for(int i=1; i < omp_get_max_threads(); ++i)
        res[0]+=res[i];

    res[0] += polynomial_cpu<vli_cpu<BaseInt, Size>, 2*Order >(gpu_product);
    return res[0];
}
#endif

template <class BaseInt, std::size_t Size, unsigned int Order>
polynomial_cpu<vli_cpu<BaseInt, Size>, 2*Order> 
inner_product_gpu( vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v1, 
               vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v2){
    std::cout<<"inner_product: single thread + CUDA"<<std::endl;
    assert(v1.size() == v2.size());
    std::size_t size_v = v1.size();
    
    polynomial_cpu<vli_cpu<BaseInt,Size>, 2*Order> res;

    std::size_t split = static_cast<std::size_t>(v1.size()*0.8);
    detail::inner_product_gpu_booster<vli_cpu<BaseInt,Size>,Order> gpu_product(v1,v2,split);

    for(std::size_t i=split ; i < size_v ; ++i){
        res += v1[i]*v2[i];
    }

    res += polynomial_cpu<vli_cpu<BaseInt, Size>, 2*Order >(gpu_product);
    
    return res;
}

} // end namespace detail
} // end namespace vli

#endif //INNER_PRODUCT_GPU_BOOSTER_HPP
