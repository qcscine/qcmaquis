#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>


namespace tools{

    boost::mt11213b rng;

    template <typename Vli>
    typename Vli::value_type rnd_digit(){
        static boost::uniform_int<typename Vli::value_type> rnd(0,vli::max_int_value<Vli>::value);
        return rnd(rng);
    }
   
    template <typename Vli>
    int rnd_valid_int(){
        static boost::uniform_int<int> rnd(0,std::abs(static_cast<int>(vli::max_int_value<Vli>::value)));
        return rnd(rng);
    }
   
    template <typename Vli>
    void vli_negate(Vli& v, int random=Vli::size){
        if(v[0]%random == 0)
            v.negate();
    }
   
    template <typename Vli>
    void fill_random(Vli& v){
        for(typename Vli::size_type i=0; i < Vli::size; ++i)
            v[i] = rnd_digit<Vli>(); 
    }
   
    template <typename Vli>
    void fill_random(Vli& v, typename Vli::size_type size){
        assert(size <= Vli::size);
        for(typename Vli::size_type i=0; i < size; ++i)
            v[i] = rnd_digit<Vli>();
    }

    template <typename Polynomial>
    void fill_poly_random(Polynomial& p){
        for(typename vli::iterator<Polynomial>::type it= p.begin(); it != p.end(); ++it){
           fill_random(*it);
           vli_negate(*it);
        }
    } 
   
    template <class Vector>
    void fill_vector_random(Vector& v){
        for(typename Vector::size_type i=0; i < v.size(); ++i)
           fill_poly_random(v[i]);
    }

   template <typename Polynomial>
   void converter(typename vli::vector_polynomial<Polynomial>& v_vli,  vli::vector_polynomial<typename vli::polynomial_multiply_type_gmp<Polynomial>::type> v_gmp ){
       typedef typename vli::polynomial_multiply_type_gmp<Polynomial>::type Polynomial_gmp;
       typedef typename vli::polynomial_multiply_type_gmp<Polynomial>::type_res Polynomial_gmp_res;
       typedef vli::vector_polynomial<Polynomial_gmp> vector_polynomial_gmp;
       typedef vli::vector_polynomial<Polynomial_gmp_res> vector_polynomial_gmp_res;

       typename Polynomial::iterator it_poly_vli;
       typename Polynomial_gmp::iterator it_poly_gmp;

       std::size_t vec_size = v_vli.size();

//       #pragma omp private(i,it_poly_vli, it_poly_gmp)
//       #pragma omp parallel for
       for(int i=0; i< vec_size; ++i){
           it_poly_vli = v_vli[i].begin();
           it_poly_gmp = v_gmp[i].begin();
           for(; it_poly_vli != v_vli[i].end(); ++it_poly_vli, ++it_poly_gmp)
               (*it_poly_gmp) = vli::detail::gmp_convert_helper<mpz_class>::apply(*it_poly_vli);
       }
   }
}
