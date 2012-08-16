
#include <benchmarks/benchmark_header.hpp>


BOOST_AUTO_TEST_CASE(TEST_NAME) // later boost mpl on the polynomial type
{
       typedef vli::polynomial_multiply_result_type<Polynomial_type_each>::type Polynomial_res;
       typedef vli::vector_polynomial<Polynomial_type_each> vector_polynomial;
       typedef vli::vector_polynomial<Polynomial_res> vector_polynomial_res;
      
       vector_polynomial v1(SIZE_VEC),v2(SIZE_VEC);
       Polynomial_res p_res;
      
       tools::fill_vector_random(v1);
       tools::fill_vector_random(v2);
      
       Timer t("CPU omp");
       t.begin();
       p_res = vli::inner_product(v1,v2);
       t.end();

} 
