
#ifndef UTILS_DMRG_RANDOM_HPP
#define UTILS_DMRG_RANDOM_HPP

#include <boost/random.hpp>

struct dmrg_random {
    typedef double value_type;
    typedef boost::mt19937 engine_t;
    typedef boost::uniform_real<value_type> uniform_dist_t;
    
    static engine_t engine;
    
    static inline value_type uniform (value_type min, value_type max) {
        uniform_dist_t dist(min, max);
        return dist(engine);
    }
    
    static inline value_type uniform () {
        return uniform(0, 1);
    }
 
};

#endif
