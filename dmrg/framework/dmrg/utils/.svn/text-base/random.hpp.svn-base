/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

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
