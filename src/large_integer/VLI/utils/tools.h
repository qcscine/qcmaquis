#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/tuple/tuple.hpp>
#include <limits>

namespace tools{
    template <class T>
    struct iterator;

    boost::mt11213b rng;

    template <typename Vli>
    typename Vli::value_type rnd_digit(){
        static boost::uniform_int<typename Vli::value_type> rnd(0,std::numeric_limits<typename Vli::value_type>::max());
        return rnd(rng);
    }
   
    template <typename Vli>
    int rnd_valid_int(){
        static boost::uniform_int<int> rnd(0,std::abs(static_cast<int>(std::numeric_limits<typename Vli::value_type>::max())));
        return rnd(rng);
    }
    // I can get a overflow during the sum of the inner product so minus 1
    template <typename Vli>
    void vli_negate(Vli& v, int random=Vli::numwords-1){
        if(v[0]%random == 0)
            v.negate();
    }
   
    template <typename Vli>
    void fill_random(Vli& v){
        for(typename Vli::size_type i=0; i < Vli::numwords-1; ++i)
            v[i] = rnd_digit<Vli>(); 
    }
   
    template <typename Vli>
    void fill_random(Vli& v, typename Vli::size_type size){
        assert(size <= Vli::numwords);
        for(typename Vli::size_type i=0; i < Vli::numwords; ++i)
            v[i] = rnd_digit<Vli>();
    }

    template <typename Polynomial>
    void fill_poly_random(Polynomial& p){
        for(typename Polynomial::iterator it= p.begin(); it != p.end(); ++it){
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
    void converter(typename vli::vector_polynomial<Polynomial>& v_vli,  vli::vector_polynomial<typename vli::polynomial_multiply_type_gmp<Polynomial>::type>& v_gmp ){
        typedef typename vli::polynomial_multiply_type_gmp<Polynomial>::type Polynomial_gmp;
        typename Polynomial::iterator it_poly_vli;
        typename Polynomial_gmp::iterator it_poly_gmp;
    
        std::size_t vec_size = v_vli.size();
    
        for(int i=0; i< vec_size; ++i){
            it_poly_vli = v_vli[i].begin();
            it_poly_gmp = v_gmp[i].begin();
            for(; it_poly_vli != v_vli[i].end(); ++it_poly_vli, ++it_poly_gmp){
                (*it_poly_gmp) = vli::detail::gmp_convert_helper<mpz_class>::apply(*it_poly_vli);
            }
        }
    }

    template <typename Polynomial>
    bool equal(typename vli::polynomial_multiply_result_type<Polynomial>::type& p_vli, typename vli::polynomial_multiply_type_gmp<Polynomial>::type_res& p_gmp){
        typedef typename vli::polynomial_multiply_type_gmp<Polynomial>::type Polynomial_gmp;
        typedef typename vli::polynomial_multiply_type_gmp<Polynomial>::type_res Polynomial_gmp_res;

        typename vli::polynomial_multiply_result_type<Polynomial>::type::iterator it_poly_vli = p_vli.begin(); 

        typename vli::polynomial_multiply_type_gmp<Polynomial>::type_res p_gmp_tmp;
        typename vli::polynomial_multiply_type_gmp<Polynomial>::type_res::iterator it_poly_gmp = p_gmp_tmp.begin(); 

        for(; it_poly_vli != p_vli.end(); ++it_poly_vli, ++it_poly_gmp)
            (*it_poly_gmp) = vli::detail::gmp_convert_helper<mpz_class>::apply(*it_poly_vli);

        return (p_gmp == p_gmp_tmp);
    }



}//end namespace tool

#endif
