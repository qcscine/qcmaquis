#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_number_traits.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "polynomial/monomial.hpp"
#include "polynomial/polynomial_gpu.hpp"
#include "polynomial/polynomial_cpu.hpp"
#include "gmpxx.h"

using vli::vli_cpu;
using vli::max_int_value;
using vli::vli_gpu;
using vli::polynomial_cpu;
using vli::monomial;

typedef boost::mpl::list<unsigned int, unsigned long int> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( constructors_test_site_monome, T, test_types)
{
    enum {poly_size = 10};
	gpu::gpu_manager* GPU;
	GPU->instance();
	
    vli::polynomial_cpu<vli_cpu<T,8>,poly_size> pa;
    vli::polynomial_cpu<vli_cpu<T,8>,poly_size> pb;
    vli::polynomial_cpu<vli_cpu<T,8>,poly_size> pc;
    
    vli::polynomial_gpu<vli_gpu<T,8>,poly_size> pagpu;
    vli::polynomial_gpu<vli_gpu<T,8>,poly_size> pcgpu;
    
    vli::polynomial_gpu<vli_gpu<T,8>,poly_size> pdgpu(pc);
    
    for(std::size_t i=0; i < poly_size; ++i){     
        for(std::size_t j=0; j < poly_size; ++j){
            for(std::size_t k=0; k < 8; ++k){
                pa(i,j)[k] = static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
                pb(i,j)[k] = static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
            }
            pagpu(i,j) = pa(i,j);
        }
    }
    
    
    vli::polynomial_gpu<vli_gpu<T,8>,poly_size> pbgpu(pb);
    
    
    /** THE ORDER IS VERY IMPORTANT, first CPU second GPU*/
    BOOST_CHECK_EQUAL(pa,pagpu);
    BOOST_CHECK_EQUAL(pb,pbgpu);
    BOOST_CHECK_EQUAL(pc,pcgpu);

	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplication_polynomial, T, test_types)
{
  	gpu::gpu_manager* GPU;
	GPU->instance();
    
    enum {poly_size = 10};

    vli::polynomial_cpu<vli_cpu<T,8>,poly_size> pa;
    vli::polynomial_cpu<vli_cpu<T,8>,poly_size> pb;
    vli::polynomial_cpu<vli_cpu<T,8>,poly_size> pc;
    
    vli::polynomial_gpu<vli_gpu<T,8>,poly_size> pagpu;
    vli::polynomial_gpu<vli_gpu<T,8>,poly_size> pbgpu;
    vli::polynomial_gpu<vli_gpu<T,8>,poly_size> pcgpu;
    
    
    for(std::size_t i=0; i < poly_size; ++i){     
        for(std::size_t j=0; j < poly_size; ++j){
            for(std::size_t k=0; k < 8; ++k){
                pa(i,j)[k] = static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
                pb(i,j)[k] = static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
            }
            pagpu(i,j) = pa(i,j);
            pbgpu(i,j) = pb(i,j);
        }
    }
    
    pc = pa*pb;
    pcgpu = pagpu*pbgpu;
        
    /** THE ORDER IS VERY IMPORTANT, first CPU second GPU*/
    BOOST_CHECK_EQUAL(pc,pcgpu);
    
	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( addition_polynomial, T, test_types)
{
    gpu::gpu_manager* GPU;
	GPU->instance();
   
    enum {poly_size = 20 };

    vli::polynomial_cpu<vli_cpu<T,8>,poly_size> pa;
    vli::polynomial_cpu<vli_cpu<T,8>,poly_size> pb;
    vli::polynomial_cpu<vli_cpu<T,8>,poly_size> pc;
    
    for(std::size_t i=0; i < poly_size; ++i){     
        for(std::size_t j=0; j < poly_size; ++j){
            for(std::size_t k=0; k < 8; ++k){
                pa(i,j)[k] = static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
                pb(i,j)[k] = static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
            }
        }
    }
    
    vli::polynomial_gpu<vli_gpu<T,8>,poly_size> pagpu(pa);
    vli::polynomial_gpu<vli_gpu<T,8>,poly_size> pbgpu(pb);
    vli::polynomial_gpu<vli_gpu<T,8>,poly_size> pcgpu(pb);

    pa+=pb;
    pagpu+=pbgpu;
    
    printf("GPU \n");
    std::cout << pagpu << std::endl;
    printf("--------------------------- \n");
    printf("CPU \n");
    std::cout << pa << std::endl;
    /** THE ORDER IS VERY IMPORTANT, first CPU second GPU*/
    BOOST_CHECK_EQUAL(pa,pagpu);
    
	GPU->instance().destructor();
    
}
