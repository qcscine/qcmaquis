/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef PARALLEL_SCHEDULER_CYCLIC_HPP
#define PARALLEL_SCHEDULER_CYCLIC_HPP

namespace parallel {

    class scheduler_cyclic {
    public:
        typedef traits::resource_iterator resource_iterator;

        resource_iterator operator()(int k) const {
            return traits::to_iterator( k % traits::size() );
        }
    };

}

#endif

