/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef PARALLEL_SCHEDULER_INPLACE_HPP
#define PARALLEL_SCHEDULER_INPLACE_HPP

namespace parallel {

    class scheduler_inplace {
    public:
        typedef traits::resource_iterator resource_iterator;

        template<typename T>
        resource_iterator operator()(const T& obj) const {
            return traits::to_iterator( traits::placement(obj) );
        }
    };

}

#endif


