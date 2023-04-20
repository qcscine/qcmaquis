/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef PARALLEL_SCHEDULER_SIZE_INDEXED_HPP
#define PARALLEL_SCHEDULER_SIZE_INDEXED_HPP

#include <algorithm>

namespace parallel {

    class scheduler_size_indexed {
    public:
        typedef traits::resource_iterator resource_iterator;

        struct index {
            typedef short size_index_type;

            template<class Matrix>
            void set(const Matrix& m){
                size_t loop_max = m.n_blocks();
                std::vector<std::pair<size_t, size_t> > tmp(loop_max);
                for(size_t k = 0; k < loop_max; ++k) 
                    tmp[k] = std::make_pair(m.left_basis()[k].second * m.right_basis()[k].second, k);
                std::sort(tmp.begin(), tmp.end());
            
                data.resize(loop_max);
                for(size_t k = 0; k < loop_max; ++k)
                    data[tmp[k].second] = k;
            }
            bool empty() const {
                return (data.size() == 0);
            }
            void insert(size_t pos, const size_index_type& val){
                if(data.size()) data.insert(data.begin() + pos, val);
            }
            void resize(size_t size){
                data.resize(size);
            }
            size_index_type& operator()(size_t k){
                if(k >= data.size()){ 
                    #ifdef USE_AMBIENT
                    printf("size_index not created!\n");
                    #endif
                    return (stub = k); 
                }
                return data[k];
            }
            friend void swap(index& x, index& y){
                std::swap(x.data, y.data);
            }
        private:
            std::vector<size_index_type> data;
            size_index_type stub;
        };

        template<class Matrix>
        scheduler_size_indexed(const Matrix& m) : i(m.size_index) {}
        resource_iterator operator()(int k) const {
            return traits::to_iterator(i(k) % traits::size());
        }
        bool propagate() const {
            return !i.empty();
        }
        index& i;
    };
}

#endif
