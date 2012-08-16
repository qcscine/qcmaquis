
#include "inner_product_header.hpp"


BOOST_AUTO_TEST_CASE_TEMPLATE( max_order_combined, T, polynomial_list_max_order_combined)
{
       typedef typename vli::polynomial_multiply_result_type<T>::type Polynomial_res;
       typedef vli::vector_polynomial<T> vector_polynomial;
       typedef vli::vector_polynomial<Polynomial_res> vector_polynomial_res;
      
       vector_polynomial v1(8),v2(8);
       Polynomial_res p1_res, p2_res;
      
       tools::fill_vector_random(v1);
       tools::fill_vector_random(v2);
      
       p1_res = vli::detail::inner_product_cpu(v1,v2);
       p2_res = vli::detail::inner_product_gpu_helper<T>::inner_product_gpu(v1,v2);
       BOOST_CHECK_EQUAL(p1_res, p2_res);
} 
