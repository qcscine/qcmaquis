/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
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
        : mpo(mpo), grid(ambient::scope::size()*s2), e2(s2, false), size1_(s1), size1(ambient::scope::size()), size2(s2) 
        {
            const std::vector<int>& except = mpo.exceptions_r;
            for(int i = 0; i < except.size(); i++) e2[except[i]] = true;
            for(int b2 = 0; b2 < size2; b2++) if(!e2[b2]) GRID(0,b2) = new block_matrix<Matrix, SymmGroup>();
        }
        void hint_left(const std::vector<block_matrix<Matrix, SymmGroup> >& t){
            for(int b1 : mpo.exceptions_l) for(int b2 = 0; b2 < size2; b2++) if(mpo.has(b1,b2)) 
                storage::hint(t[b1], ambient::scope::permute(b2,mpo.placement_r));
        }
        void hint_right(Boundary<Matrix, SymmGroup> const & t){
            for(int b1 = 0; b1 < size1_; b1++)
            for(int b2 = 0; b2 < size2; b2++){
                if(!mpo.has(b1,b2)) continue;
                storage::hint(t[b2], where(b1,b2));
            }
        }
        ambient::scope::const_iterator where(size_t b1, size_t b2){
            if(e2[b2]) return ambient::scope::permute(b1,mpo.placement_l);
            return ambient::scope::permute(b2,mpo.placement_r);
        }
        block_matrix<Matrix, SymmGroup>& operator()(size_t b1, size_t b2){
            if(e2[b2]){
                block_matrix<Matrix, SymmGroup>*& el = GRID(ambient::which(),b2);
                if(el == NULL) el = new block_matrix<Matrix, SymmGroup>();
                return *el;
            }
            return *GRID(0,b2);
        }
        void multiply_column(size_t b2, const block_matrix<Matrix, SymmGroup>& rhs){
            for(int b1 = 0; b1 < size1; b1++){
                if(GRID(b1,b2) == NULL) continue;
                select_scope(e2[b2] ? (ambient::scope::begin()+b1) : where(b1,b2));
                block_matrix<Matrix, SymmGroup> res;
                gemm(*GRID(b1,b2), rhs, res);
                swap(*GRID(b1,b2), res);
                if(!e2[b2]) return;
            }
        }
        void multiply_column_trans(size_t b2, const block_matrix<Matrix, SymmGroup>& rhs){
            block_matrix<Matrix, SymmGroup> tmp;
            block_matrix<Matrix, SymmGroup> red = reduce_column(b2);
            select_scope(ambient::scope::permute(b2,mpo.placement_r));
            gemm(transpose(red), rhs, tmp);
            swap(*GRID(0,b2), tmp);
        }
        Boundary<Matrix, SymmGroup> make_boundary(){
            Boundary<Matrix, SymmGroup> ret; ret.resize(size2);
            for(int b2 = 0; b2 < size2; b2++) swap(ret[b2], reduce_column(b2));
            return ret;
        }
        block_matrix<Matrix, SymmGroup>& reduce_column(size_t b2){
            if(!e2[b2]) return *GRID(0,b2);
            ambient::scope::const_iterator r = ambient::scope::permute(b2,mpo.placement_r);

            std::vector< std::pair< block_matrix<Matrix, SymmGroup>*, ambient::scope::const_iterator > > rvector;
            for(int b1 = *r; b1 < size1; b1++) if(GRID(b1,b2) != NULL) rvector.push_back(std::make_pair(GRID(b1,b2), ambient::scope::begin()+b1));
            for(int b1 = 0; b1 < *r; b1++)     if(GRID(b1,b2) != NULL) rvector.push_back(std::make_pair(GRID(b1,b2), ambient::scope::begin()+b1));
            for(int b1 = 0; b1 < size1; b1++) GRID(b1,b2) = NULL;

            GRID(0,b2) = ambient::reduce(rvector, [](std::pair<block_matrix<Matrix, SymmGroup>*, ambient::scope::const_iterator>& dst_pair, 
                                                     std::pair<block_matrix<Matrix, SymmGroup>*, ambient::scope::const_iterator>& src_pair){
                                                         block_matrix<Matrix, SymmGroup>& dst = *dst_pair.first;
                                                         block_matrix<Matrix, SymmGroup>& src = *src_pair.first;
                                                         select_scope(dst_pair.second);
                                                         for(size_t k = 0; k < src.n_blocks(); ++k)
                                                         dst.match_and_add_block(src[k],
                                                                                 src.left_basis()[k].first, 
                                                                                 src.right_basis()[k].first);
                                                     }).first;
            e2[b2] = false;
            if(rvector[0].second != r) storage::migrate(*rvector[0].first, r);
            for(int i = 1; i < rvector.size(); i++) delete rvector[i].first;
            return *GRID(0,b2);
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
        block_matrix<Matrix, SymmGroup>& reduce(){
            std::vector< std::pair< block_matrix<Matrix, SymmGroup>*, ambient::scope::const_iterator > > rvector;
            std::vector< std::pair< block_matrix<Matrix, SymmGroup>*, ambient::scope::const_iterator >* > rvector_global;
        
            for(int b2 = 0; b2 < size2; b2++)
            for(int b1 = 0; b1 < size1; b1++){
                if(GRID(b1,b2) == NULL) continue;
                ambient::scope::const_iterator owner = e2[b2] ? (ambient::scope::begin()+b1) : ambient::scope::permute(b2,mpo.placement_r);
                rvector.push_back(std::make_pair(GRID(b1,b2), owner));
                GRID(b1,b2) = NULL;
            }
        
            std::sort(rvector.begin(), rvector.end(), [](const std::pair<block_matrix<Matrix, SymmGroup>*, ambient::scope::const_iterator>& a, 
                                                         const std::pair<block_matrix<Matrix, SymmGroup>*, ambient::scope::const_iterator>& b){ 
                                                             return a.second < b.second; 
                                                         });
            int i = 0;
            while(i < rvector.size()){
                ambient::scope::const_iterator owner = rvector[i].second;
                rvector_global.push_back(&rvector[i++]);
                select_scope(owner);
                while(rvector[i].second == owner && i < rvector.size()){
                    for(size_t k = 0; k < rvector[i].first->n_blocks(); ++k)
                    rvector_global.back()->first->match_and_add_block((*rvector[i].first)[k], 
                                                                      rvector[i].first->left_basis()[k].first, 
                                                                      rvector[i].first->right_basis()[k].first);
                    delete rvector[i].first;
                    i++;
                }
            }

            std::vector<Matrix*> blocks;
            std::vector<typename SymmGroup::charge> c1;
            std::vector<typename SymmGroup::charge> c2;
            std::vector<ambient::scope::const_iterator> owners;
            block_matrix<Matrix, SymmGroup>* skeleton = new block_matrix<Matrix, SymmGroup>();
            block_matrix<Matrix, SymmGroup>* res = new block_matrix<Matrix, SymmGroup>();

            for(size_t n = 0; n < rvector_global.size(); ++n){
                block_matrix<Matrix, SymmGroup>& src = *rvector_global[n]->first;
                for(size_t k = 0; k < src.n_blocks(); ++k){
                    if(!skeleton->has_block(src.left_basis()[k].first, src.right_basis()[k].first)){
                        skeleton->insert_block(new Matrix(), src.left_basis()[k].first, src.right_basis()[k].first);
                    }
                    blocks.push_back(&src[k]);
                    c1.push_back(src.left_basis()[k].first);
                    c2.push_back(src.right_basis()[k].first);
                    owners.push_back(rvector_global[n]->second);
                }
            }
            std::vector< std::vector<std::pair<Matrix*,ambient::scope::const_iterator> > > rblocks;
            size_t max_stride = 0;
            for(size_t k = 0; k < skeleton->n_blocks(); ++k){
                auto tc1 = skeleton->left_basis()[k].first; 
                auto tc2 = skeleton->right_basis()[k].first;
                std::vector<std::pair<Matrix*,ambient::scope::const_iterator> > rblocks_part;
                for(size_t n = 0; n < blocks.size(); n++){
                    if(tc1 == c1[n] && tc2 == c2[n]) rblocks_part.push_back(std::make_pair(blocks[n], owners[n]));
                }
                ambient::scope::const_iterator root = ambient::scope::balance(k,skeleton->n_blocks());
                std::sort(rblocks_part.begin(), rblocks_part.end(), [root](const std::pair<Matrix*, ambient::scope::const_iterator>& a,
                                                                           const std::pair<Matrix*, ambient::scope::const_iterator>& b){
                                                                               return (ambient::num_procs() + *a.second - *root) % ambient::num_procs()
                                                                                    < (ambient::num_procs() + *b.second - *root) % ambient::num_procs();
                                                                           });
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
                        select_scope(dst_pair.second);
                        Matrix& src = *src_pair.first;
                        Matrix& dst = *dst_pair.first;
                        
                        if(num_rows(src) == num_rows(dst) && num_cols(src) == num_cols(dst)){
                        }else if(num_rows(src) > num_rows(dst) && num_cols(src) > num_cols(dst))
                            resize(dst, num_rows(src), num_cols(src));
                        else{
                            size_t maxrows = std::max(num_rows(src), num_rows(dst));
                            size_t maxcols = std::max(num_cols(src), num_cols(dst));
                            resize(dst, maxrows, maxcols);
                            resize(src, maxrows, maxcols);
                        }
                        dst += src;
                        Matrix tmp; src.swap(tmp);
                    }
                }
                ambient::sync();
            }
            for(size_t k = 0; k < skeleton->n_blocks(); ++k){
                auto tc1 = skeleton->left_basis()[k].first; 
                auto tc2 = skeleton->right_basis()[k].first;
                res->insert_block(*rblocks[k][0].first, tc1, tc2);
            }
            for(size_t n = 0; n < rvector_global.size(); ++n) delete rvector_global[n]->first;
            delete skeleton;
            GRID(0,0) = res;
            return *res;
        }
        MPOTensor<Matrix, SymmGroup> const & mpo;
        mutable std::vector< block_matrix<Matrix, SymmGroup>* > grid;
        std::vector<bool> e2;
        size_t size1_; // actual
        size_t size1;  // optimized
        size_t size2;
        #undef GRID
    };
}

#endif
