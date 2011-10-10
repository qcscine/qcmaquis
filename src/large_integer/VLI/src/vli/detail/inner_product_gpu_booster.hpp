#ifndef INNER_PRODUCT_GPU_BOOSTER_HPP
#define INNER_PRODUCT_GPU_BOOSTER_HPP

#include "vli/utils/gpu_manager.h"
//#include "utils/timings.h"

#include "vli/detail/kernels_gpu.h"
#include <iostream>

namespace vli
{
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
    template <typename VectorPolynomial>
    class inner_product_gpu_booster
    {
      private:
        typedef typename VectorPolynomial::value_type   polynomial_type;
        typedef typename polynomial_type::value_type    vli_type;
        typedef typename vli_type::value_type           base_int_type;
        enum {element_size = polynomial_type::max_order * polynomial_type::max_order * vli_type::size };
        enum {num_threads = 256}; // TODO change
      public:
        typedef std::size_t size_type;

        explicit inner_product_gpu_booster(VectorPolynomial const& v1, VectorPolynomial const& v2, size_type partsize)
        : v1_(partsize*element_size), v2_(partsize*element_size), tmp_(partsize*element_size)
        {
            assert(partsize <= v1.size());
            assert(partsize <= v2.size());
            gpu::cu_check_error(cudaMemcpyAsync((void*)v1_.p(),(void*)&v1[0],partsize*element_size*sizeof(base_int_type),cudaMemcpyHostToDevice), __LINE__);
            gpu::cu_check_error(cudaMemcpyAsync((void*)v2_.p(),(void*)&v2[0],partsize*element_size*sizeof(base_int_type),cudaMemcpyHostToDevice),__LINE__);
            gpu::cu_check_error(cudaMemset((void*)tmp_.p(),0,partsize*element_size*sizeof(base_int_type)),__LINE__);
            

            // run/queue the inner product computation
//            TimerCuda A("produit");
//            A.begin();
            inner_product_vector(vli_size_tag<vli_type::size>(),polynomial_type::max_order, partsize,v1_.p(), v2_.p(), tmp_.p(),num_threads);
            gpu::cu_check_error(cudaGetLastError(),__LINE__);
            vector_reduction_inplace(vli_size_tag<vli_type::size>(), polynomial_type::max_order, partsize, tmp_.p());
            gpu::cu_check_error(cudaGetLastError(),__LINE__);
        }
        
        operator polynomial_type() const
        {
            polynomial_type poly;
            gpu::cu_check_error(cudaMemcpy((void*)&poly(0,0),(void*)tmp_.p(),element_size*sizeof(base_int_type),cudaMemcpyDeviceToHost),__LINE__);
            return poly;
        }

      private:
        gpu_memblock<base_int_type> v1_;
        gpu_memblock<base_int_type> v2_;
        gpu_memblock<base_int_type> tmp_;
    };
}

template <class BaseInt, std::size_t Size>
class vli_cpu;

template <class Vli, unsigned int Order>
class polynomial_cpu;

template <class Polynomial>
class vector_polynomial_cpu;

namespace detail
{
#ifdef _OPENMP
template <class BaseInt, std::size_t Size, unsigned int Order>
polynomial_cpu<vli_cpu<BaseInt, Size>, Order> 
inner_product_openmp_gpu( vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v1, 
               vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v2){

//    std::cout<<"inner_product: OpenMP +CUDA"<<std::endl;
    assert(v1.size() == v2.size());
    std::size_t size_v = v1.size();
    polynomial_cpu<vli_cpu<BaseInt, Size>, Order>  res[omp_get_max_threads()];
    
    std::size_t split = static_cast<std::size_t>(v1.size()*0.1);
    detail::inner_product_gpu_booster<vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> > > gpu_product(v1,v2,split);

    #pragma omp parallel for
    for(std::size_t i=split ; i < size_v ; ++i){
        res[omp_get_thread_num()] += v1[i]*v2[i];
    }

    for(int i=1; i < omp_get_max_threads(); ++i)
        res[0]+=res[i];

    res[0] += polynomial_cpu<vli_cpu<BaseInt, Size>, Order >(gpu_product);
    
    return res[0];
}
#endif
template <class BaseInt, std::size_t Size, unsigned int Order>
polynomial_cpu<vli_cpu<BaseInt, Size>, Order> 
inner_product_gpu( vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v1, 
               vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> >  const& v2){
//    std::cout<<"inner_product: single thread + CUDA"<<std::endl;
    assert(v1.size() == v2.size());
    std::size_t size_v = v1.size();
    
    polynomial_cpu<vli_cpu<BaseInt,Size>, Order> res;

    std::size_t split = static_cast<std::size_t>(v1.size()*0.1);
    detail::inner_product_gpu_booster<vector_polynomial_cpu<polynomial_cpu<vli_cpu<BaseInt, Size>, Order> > > gpu_product(v1,v2,split);

    for(std::size_t i=split ; i < size_v ; ++i){
        res += v1[i]*v2[i];
    }

    res += polynomial_cpu<vli_cpu<BaseInt, Size>, Order >(gpu_product);
    
    return res;
}

} // end namespace detail
} // end namespace vli

#endif //INNER_PRODUCT_GPU_BOOSTER_HPP
