/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef PARALLEL_SCHEDULER_BALANCED_HPP
#define PARALLEL_SCHEDULER_BALANCED_HPP

namespace parallel {

    class scheduler_balanced {
    public:
        typedef traits::resource_iterator resource_iterator;

        template<class Matrix>
        scheduler_balanced(const Matrix& m) : max_k(m.n_blocks()) {}
        scheduler_balanced(size_t max) : max_k(max) {}
        resource_iterator operator()(int k) const {
            return traits::balance(k,max_k);
        }
        bool propagate() const {
            return false;
        }
    protected:
        int max_k;
    };

}

#endif
