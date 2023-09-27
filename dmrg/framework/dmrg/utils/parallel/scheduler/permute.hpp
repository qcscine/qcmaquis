/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef PARALLEL_SCHEDULER_PERMUTE_HPP
#define PARALLEL_SCHEDULER_PERMUTE_HPP

namespace parallel {

    class scheduler_permute {
    public:
        typedef traits::resource_iterator resource_iterator;

        scheduler_permute(const std::vector<int>& s, int gran = 1) : permutation(s), granularity(gran) {}

        resource_iterator operator()(int b) const {
            return traits::permute(b, permutation, granularity);
        }
    private:
        const std::vector<int>& permutation;
        int granularity;
    };

}

#endif


