/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2012 by Alexandr Kosenkov <alex.kosenkov@gmail.com>
 *                            Timothee Ewart <timothee.ewart@gmail.com>
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

#ifndef CONTRACTIONS_IMPL_AMBIENT_HPP
#define CONTRACTIONS_IMPL_AMBIENT_HPP

namespace contraction {

    template<class Matrix, class SymmGroup>
    class ContractionGrid {
    public:
        #define GRID(b1,b2) grid[b1+b2*size1]
        static void check_exceptions(MPOTensor<Matrix, SymmGroup> const & mpo, int s1, int s2){
            std::vector<int> ex;
            for(int b2 = 0; b2 < s2; b2++){
                bool present = false;
                for(int b1 = 0; b1 < s1; b1++) if(mpo.has(b1,b2)){
                    if(present){ ex.push_back(b2); break; }
                    present = true;
                }
            }
            if(ex == mpo.exceptions_r) return;
            int size = ex.size();
            if(ex.size() != mpo.exceptions_r.size()){
                std::cout << "Error: sizes are different!\n";
                size = std::min(ex.size(), mpo.exceptions_r.size());
            }
            for(int i = 0; i < size; i++) if(ex[i] != mpo.exceptions_r[i])
                std::cout << "Error: exceptions are different!\n";
        }
       ~ContractionGrid(){
            for(int i = 0; i < size1*size2; i++) delete grid[i];
        }
        ContractionGrid(MPOTensor<Matrix, SymmGroup> const & mpo, size_t s1, size_t s2)
        : mpo(mpo), 
          size1_(s1),
          size2(s2),
          e2(s2, false)
        {
            granularity = parallel::groups_granularity;
            if(granularity > ambient::scope::size()) granularity = 1;
            size1 = ambient::scope::size()/granularity;
            grid.resize(size1*size2);
            const std::vector<int>& except = mpo.exceptions_r;
            for(int i = 0; i < except.size(); i++) e2[except[i]] = true;
            for(int b2 = 0; b2 < size2; b2++) if(!e2[b2]) GRID(0,b2) = new block_matrix<Matrix, SymmGroup>();
        }
        static void iterate_reduction_layout(int start, int end){
            reduction_iteration = start;
            reduction_max_iteration = end;
        }
        void index_sizes(size_t b2){
            for(int b1 = 0; b1 < size1; b1++){
                if(GRID(b1,b2) == NULL) continue;
                GRID(b1,b2)->index_sizes();
            }
        }
        ambient::scope::const_iterator where(size_t b1, size_t b2){
            if(e2[b2]) return ambient::scope::permute(b1,mpo.placement_l,granularity);
            return ambient::scope::permute(b2,mpo.placement_r,granularity);
        }
        block_matrix<Matrix, SymmGroup>& operator()(size_t b1, size_t b2){
            if(e2[b2]){
                block_matrix<Matrix, SymmGroup>*& el = GRID((where(b1,b2)-ambient::scope::begin())/granularity,b2);
                if(el == NULL) el = new block_matrix<Matrix, SymmGroup>();
                return *el;
            }
            return *GRID(0,b2);
        }
        void multiply_column(size_t b2, const block_matrix<Matrix, SymmGroup>& rhs){
            for(int b1 = 0; b1 < size1; b1++){
                if(GRID(b1,b2) == NULL) continue;
                auto it = e2[b2] ? (ambient::scope::begin()+b1*granularity) : where(b1,b2);
                parallel::guard group(it, granularity);
                block_matrix<Matrix, SymmGroup> res;
                gemm(*GRID(b1,b2), rhs, res, parallel::scheduler_size_indexed(*GRID(b1,b2)));
                swap(*GRID(b1,b2), res);
                if(!e2[b2]) return;
            }
        }
        void multiply_column_trans(size_t b2, const block_matrix<Matrix, SymmGroup>& rhs){
            block_matrix<Matrix, SymmGroup> tmp;
            auto red = transpose(reduce_column(b2));
            parallel::guard group(ambient::scope::permute(b2, mpo.placement_r, granularity), granularity);
            gemm(red, rhs, tmp, parallel::scheduler_size_indexed(red));
            swap(*GRID(0,b2), tmp);
        }
        Boundary<Matrix, SymmGroup> make_boundary(){
            Boundary<Matrix, SymmGroup> ret; ret.resize(size2);
            for(int b2 = 0; b2 < size2; b2++) swap(ret[b2], reduce_column(b2));
            return ret;
        }
        block_matrix<Matrix, SymmGroup>& reduce_column(size_t b2){
            if(!e2[b2]) return *GRID(0,b2);
            ambient::scope::const_iterator center = ambient::scope::permute(b2,mpo.placement_r,granularity);
            int r = (center-ambient::scope::begin())/granularity;

            std::vector< std::pair< block_matrix<Matrix, SymmGroup>*, ambient::scope::const_iterator > > rvector;
            for(int b1 = r; b1 < size1; b1++) if(GRID(b1,b2) != NULL) rvector.push_back(std::make_pair(GRID(b1,b2), ambient::scope::begin()+b1*granularity));
            for(int b1 = 0; b1 < r; b1++)     if(GRID(b1,b2) != NULL) rvector.push_back(std::make_pair(GRID(b1,b2), ambient::scope::begin()+b1*granularity));
            for(int b1 = 0; b1 < size1; b1++) GRID(b1,b2) = NULL;

            GRID(0,b2) = ambient::reduce(rvector, [this](std::pair<block_matrix<Matrix, SymmGroup>*, ambient::scope::const_iterator>& dst_pair, 
                                                         std::pair<block_matrix<Matrix, SymmGroup>*, ambient::scope::const_iterator>& src_pair){
                                                             parallel::scheduler_cyclic scheduler;

                                                             block_matrix<Matrix, SymmGroup>& dst = *dst_pair.first;
                                                             block_matrix<Matrix, SymmGroup>& src = *src_pair.first;
                                                             parallel::guard group(dst_pair.second, granularity);
                                                             for(size_t k = 0; k < src.n_blocks(); ++k){
                                                                 size_t index = dst.find_block(src.basis().left_charge(k), src.basis().right_charge(k));
                                                                 // alternative: index = src.size_index(k);
                                                                 parallel::guard proc(scheduler(index));
                                                                 dst.match_and_add_block(src[k],
                                                                                         src.basis().left_charge(k), 
                                                                                         src.basis().right_charge(k));
                                                             }
                                                         }).first;
            e2[b2] = false;

            rvector[0].first->index_sizes();
            parallel::guard group(center, granularity);
            storage::migrate(*rvector[0].first, parallel::scheduler_size_indexed(*rvector[0].first));

            for(int i = 1; i < rvector.size(); i++) delete rvector[i].first;
            return *GRID(0,b2);
        }
        block_matrix<Matrix, SymmGroup>& reduce(){
            std::vector< std::pair< block_matrix<Matrix, SymmGroup>*, ambient::scope::const_iterator > > rvector;
        
            for(int b2 = 0; b2 < size2; b2++)
            for(int b1 = 0; b1 < size1; b1++){
                if(GRID(b1,b2) == NULL) continue;
                ambient::scope::const_iterator owner = e2[b2] ? (ambient::scope::begin()+b1*granularity) : ambient::scope::permute(b2,mpo.placement_r,granularity);
                rvector.push_back(std::make_pair(GRID(b1,b2), owner));
                GRID(b1,b2) = NULL;
            }
        
            std::sort(rvector.begin(), rvector.end(), [](const std::pair<block_matrix<Matrix, SymmGroup>*, ambient::scope::const_iterator>& a, 
                                                         const std::pair<block_matrix<Matrix, SymmGroup>*, ambient::scope::const_iterator>& b){ 
                                                             return a.second < b.second; 
                                                         });
            std::vector<Matrix*> blocks;
            std::vector<typename SymmGroup::charge> c1;
            std::vector<typename SymmGroup::charge> c2;
            std::vector<ambient::scope::const_iterator> owners;
            block_matrix<Matrix, SymmGroup>* skeleton = new block_matrix<Matrix, SymmGroup>();
            block_matrix<Matrix, SymmGroup>* res = new block_matrix<Matrix, SymmGroup>();

            for(size_t n = 0; n < rvector.size(); ++n){
                block_matrix<Matrix, SymmGroup>& src = *rvector[n].first;
                for(size_t k = 0; k < src.n_blocks(); ++k){
                    if(!skeleton->has_block(src.basis().left_charge(k), src.basis().right_charge(k))){
                        skeleton->insert_block(new Matrix(), src.basis().left_charge(k), src.basis().right_charge(k));
                    }
                    blocks.push_back(&src[k]);
                    c1.push_back(src.basis().left_charge(k));
                    c2.push_back(src.basis().right_charge(k));
                    owners.push_back(rvector[n].second + src.size_index(k)%granularity);
                }
            }

            res->index_iter(reduction_iteration, reduction_max_iteration);
            skeleton->index_iter(reduction_iteration, reduction_max_iteration);
            parallel::scheduler_balanced_iterative scheduler(*skeleton);
            ++reduction_iteration %= reduction_max_iteration;

            std::vector< std::vector<std::pair<Matrix*,ambient::scope::const_iterator> > > rblocks;
            size_t max_stride = 0;
            for(size_t k = 0; k < skeleton->n_blocks(); ++k){
                auto tc1 = skeleton->basis().left_charge(k); 
                auto tc2 = skeleton->basis().right_charge(k);
                std::vector<std::pair<Matrix*,ambient::scope::const_iterator> > rblocks_part;
                for(size_t n = 0; n < blocks.size(); n++){
                    if(tc1 == c1[n] && tc2 == c2[n]) rblocks_part.push_back(std::make_pair(blocks[n], owners[n]));
                }
                ambient::scope::const_iterator root = scheduler(k);
                std::sort(rblocks_part.begin(), rblocks_part.end(), [root](const std::pair<Matrix*, ambient::scope::const_iterator>& a,
                                                                           const std::pair<Matrix*, ambient::scope::const_iterator>& b){
                                                                               return (ambient::num_procs() + *a.second - *root) % ambient::num_procs()
                                                                                    < (ambient::num_procs() + *b.second - *root) % ambient::num_procs();
                                                                           });
                rblocks_part[0].second = root;
                rblocks.push_back(rblocks_part);
                if(rblocks_part.size() > max_stride) 
                    max_stride = rblocks_part.size();
            }
            for(int stride = 1; stride < max_stride; stride *= 2){
                for(size_t n = 0; n < rblocks.size(); ++n){
                    auto& rblocks_part = rblocks[n];
                    if(stride < rblocks_part.size())
                    for(int k = stride; k < rblocks_part.size(); k += stride*2){
                        std::pair< Matrix*, ambient::scope::const_iterator >& dst_pair = rblocks_part[k-stride];
                        std::pair< Matrix*, ambient::scope::const_iterator >& src_pair = rblocks_part[k];
                        parallel::guard proc(dst_pair.second);
                        Matrix& src = *src_pair.first;
                        Matrix& dst = *dst_pair.first;
                        
                        if(num_rows(src) > num_rows(dst) && num_cols(src) > num_cols(dst))
                            resize(dst, num_rows(src), num_cols(src));
                        else if(num_rows(src) != num_rows(dst) || num_cols(src) != num_cols(dst)){
                            size_t maxrows = std::max(num_rows(src), num_rows(dst));
                            size_t maxcols = std::max(num_cols(src), num_cols(dst));
                            resize(dst, maxrows, maxcols);
                            resize(src, maxrows, maxcols);
                        }
                        dst += src;
                        Matrix tmp; src.swap(tmp);
                    }
                }
                parallel::sync();
            }
            for(size_t k = 0; k < skeleton->n_blocks(); ++k){
                auto tc1 = skeleton->basis().left_charge(k); 
                auto tc2 = skeleton->basis().right_charge(k);
                res->insert_block(*rblocks[k][0].first, tc1, tc2);
            }
            for(size_t n = 0; n < rvector.size(); ++n) delete rvector[n].first;
            delete skeleton;
            GRID(0,0) = res;
            return *res;
        }
        void print_distribution(){
            if(!ambient::master()) return;
            double total = 0;
            for(int b2 = 0; b2 < size2; b2++)
            for(int b1 = 0; b1 < size1; b1++) if(GRID(b1,b2) != NULL) 
                total += GRID(b1,b2)->num_elements();
            printf("%.2f GB:", total*sizeof(typename Matrix::value_type)/1024/1024/1024);
            for(int p = 0; p < ambient::num_procs(); ++p){
                double part = 0;
                for(int b2 = 0; b2 < size2; b2++)
                for(int b1 = 0; b1 < size1; b1++) if(GRID(b1,b2) != NULL)
                for(int i = 0; i < (*GRID(b1,b2)).n_blocks(); ++i){
                    if(!ambient::weak((*GRID(b1,b2))[i][0]) && ambient::get_owner((*GRID(b1,b2))[i][0]) == p)
                        part += num_rows((*GRID(b1,b2))[i])*num_cols((*GRID(b1,b2))[i]);
                }
                printf(" %.1f%%", 100*part/total);
            }
            printf("\n");
        }
        MPOTensor<Matrix, SymmGroup> const & mpo;
        mutable std::vector< block_matrix<Matrix, SymmGroup>* > grid;
        static int reduction_max_iteration;
        static int reduction_iteration;
        std::vector<bool> e2;
        size_t granularity;
        size_t size1_; // actual
        size_t size1;  // optimized
        size_t size2;
        #undef GRID
    };

    template<class Matrix, class SymmGroup> int ContractionGrid<Matrix, SymmGroup>::reduction_iteration = 0;
    template<class Matrix, class SymmGroup> int ContractionGrid<Matrix, SymmGroup>::reduction_max_iteration = 1;
}

#endif
