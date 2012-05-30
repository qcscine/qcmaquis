
#include <benchmarks/benchmark_header.hpp>

using namespace vli::bench;

BOOST_AUTO_TEST_CASE(TEST_NAME) // later boost mpl on the polynomial type
{
//    assert(boost::unit_test::framework::master_test_suite().argc == 1);
//    int SizeVector = atoi(boost::unit_test::framework::master_test_suite().argv[0]); 
    int SizeVector = 512;
     
    vector_type v1(SizeVector);
    vector_type v2(SizeVector);

    polynomial_type_double result_cpu  ;
       
    Timer ti("Timer init");
    ti.begin();
    fill_vector_random(v1,VLI_SIZE);
    fill_vector_random(v2,VLI_SIZE-1);

    fill_vector_negate(v1,VLI_SIZE);
    fill_vector_negate(v2,VLI_SIZE-1);
    ti.end();
   
    std::string name("InnerProduct");
   
    name += "_VLI_SIZE_" + boost::lexical_cast<std::string>(VLI_SIZE);
    name += "_POLY_ORDER_" + boost::lexical_cast<std::string>(ORDER_POLY); // ascii art ~~' 

    Timer t(name);
    t.begin();
      result_cpu = vli::detail::inner_product_plain(v1,v2);
    t.end();
} 
