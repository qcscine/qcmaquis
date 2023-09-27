/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef PARALLEL_SCHEDULER_NOP_HPP
#define PARALLEL_SCHEDULER_NOP_HPP

namespace parallel {

    class scheduler_nop {
    public:
        typedef traits::resource_iterator_nop resource_iterator;

        resource_iterator operator()(int k) const {
            return resource_iterator();
        }
        bool propagate() const {
            return false;
        }
    };

}

#endif
