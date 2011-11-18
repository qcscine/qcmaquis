
#ifndef UTILS_DMRG_RANDOM_HPP
#define UTILS_DMRG_RANDOM_HPP

#include <boost/random.hpp>

struct dmrg_random {
    typedef double value_type;
    typedef boost::mt19937 engine_t;
    typedef boost::uniform_real<value_type> uniform_dist_t;
    typedef boost::normal_distribution<value_type> normal_dist_t;
    typedef boost::poisson_distribution<value_type> poisson_dist_t;
    
    static engine_t engine;

    // Uniform distribution
    static inline value_type uniform (value_type min, value_type max) {
        uniform_dist_t dist(min, max);
        return dist(engine);
    }
    
    static inline value_type uniform () {
        return uniform(0, 1);
    }

    
    // Normal distribution
    static inline value_type normal (value_type mean, value_type sigma) {
        normal_dist_t dist(mean, sigma);
        return dist(engine);
    }
    
    static inline value_type normal () {
        return normal(0, 1);
    }

    
    // Poisson distribution
    /*
    static inline value_type poisson (value_type mean) {
        poisson_dist_t dist(mean);
        return dist(engine);
    }
    
    static inline value_type poisson () {
        return poisson(1);
    }
*/
    
};

#endif
