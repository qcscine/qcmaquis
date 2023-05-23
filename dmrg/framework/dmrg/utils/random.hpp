/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_DMRG_RANDOM_HPP
#define UTILS_DMRG_RANDOM_HPP

#include <random>

struct dmrg_random {
    using value_type = double;
    using engine_t = std::mt19937;
    using uniform_dist_t = std::uniform_real_distribution<>;
    using normal_dist_t = std::normal_distribution<>;
    using poisson_dist_t = std::poisson_distribution<value_type>;
    
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
