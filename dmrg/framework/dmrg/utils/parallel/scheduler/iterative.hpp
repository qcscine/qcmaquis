/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Alexandr Kosenkov <alex.kosenkov@gmail.com>
 *                         by Michele Dolfi <dolfim@phys.ethz.ch>
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
