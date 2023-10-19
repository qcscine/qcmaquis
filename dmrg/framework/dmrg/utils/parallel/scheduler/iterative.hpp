/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef PARALLEL_SCHEDULER_ITERATIVE_HPP
#define PARALLEL_SCHEDULER_ITERATIVE_HPP

namespace parallel {

    class scheduler_balanced_iterative : public scheduler_balanced {
    public:
        typedef traits::resource_iterator resource_iterator;

        struct index {
            index() : max_iterations(1), iteration(0) {}
            void set(int i, int max){
                iteration = i;
                max_iterations = max;
            }
            friend void swap(index& x, index& y){
                std::swap(x.max_iterations, y.max_iterations);
                std::swap(x.iteration,      y.iteration);
            }
            int iteration;
            int max_iterations;
        };

        template<class Matrix>
        scheduler_balanced_iterative(const Matrix& m) : scheduler_balanced(m) {
            iteration  = m.iter_index.iteration;
            chunks_num = m.iter_index.max_iterations;
            chunk_size = std::max(1, this->max_k/chunks_num);
            group_size = std::max((size_t)1, traits::size()/chunks_num);
        }

        resource_iterator operator()(int k) const {
            if(chunks_num > 1){
                int group_id = (chunks_num + (int)(k/chunk_size) - iteration) % chunks_num;
                return traits::to_iterator((group_id*group_size + k % group_size) % traits::size());
            }
            return scheduler_balanced::operator()(k);
        }
    private:
        int iteration;
        int chunks_num;
        int chunk_size;
        int group_size;
    };

}

#endif
