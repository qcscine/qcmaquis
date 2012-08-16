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
    void vli_negate(Vli& v, int random=Vli::size-1){
        if(v[0]%random == 0)
            v.negate();
    }
   
    template <typename Vli>
    void fill_random(Vli& v){
        for(typename Vli::size_type i=0; i < Vli::size-1; ++i)
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
}
