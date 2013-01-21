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
        return   rnd(rng);
    }
   
    template <typename Vli>
    int rnd_valid_int(){
        static boost::uniform_int<int> rnd(0,std::abs(static_cast<int>(std::numeric_limits<typename Vli::value_type>::max())));
        return  rnd(rng);
    }
    // I can get a overflow during the sum of the inner product so minus 1
    template <typename Vli>
    void vli_negate(Vli& v, int random=Vli::numwords-1){
        if(v[0]%random == 0)
            v.negate();
    }
   
    template <typename Vli>
    void fill_random(Vli& v){
        for(int i=0; i < Vli::numwords-1; ++i)
            v[i] = rnd_digit<Vli>();
    }
   
    template <typename Vli>
    void fill_random(Vli& v, typename Vli::size_type size){
        assert(size <= Vli::numwords);
        for(int i=0; i < size; ++i)
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

}//end namespace tool

#endif
