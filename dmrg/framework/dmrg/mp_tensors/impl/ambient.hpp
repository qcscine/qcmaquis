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

#ifndef MAQUIS_MP_TENSORS_IMPL_AMBIENT_HPP
#define MAQUIS_MP_TENSORS_IMPL_AMBIENT_HPP

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
        : mpo(mpo), grid(ambient::num_workers()*s2), e2(s2, false), size1(ambient::num_workers()), size2(s2) 
        {
            const std::vector<int>& except = mpo.exceptions_r;
            for(int i = 0; i < except.size(); i++) e2[except[i]] = true;
            for(int b2 = 0; b2 < size2; b2++) if(!e2[b2]) GRID(0,b2) = new block_matrix<Matrix, SymmGroup>();
        }
        int where(size_t b1, size_t b2){
            if(e2[b2]) return ambient::scope::permute(b1,mpo.placement_l);
            return ambient::scope::permute(b2,mpo.placement_r);
        }
        block_matrix<Matrix, SymmGroup>& operator()(size_t b1, size_t b2){
            ambient::get_scope().set(where(b1,b2));
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
                ambient::get_scope().set(e2[b2] ? b1 : where(b1,b2));
                block_matrix<Matrix, SymmGroup> res;
                gemm(*GRID(b1,b2), rhs, res);
                *GRID(b1,b2) = res;
                if(!e2[b2]) return;
            }
        }
        void multiply_column_trans(size_t b2, const block_matrix<Matrix, SymmGroup>& rhs){
            block_matrix<Matrix, SymmGroup> tmp;
            gemm(transpose(reduce_column(b2)), rhs, tmp);
            *GRID(0,b2) = tmp;
        }
        Boundary<Matrix, SymmGroup> make_boundary(){
            select_proc(0);
            Boundary<Matrix, SymmGroup> ret; ret.resize(size2);
            for(int b2 = 0; b2 < size2; b2++) ret[b2] = reduce_column(b2);
            return ret;
        }
        block_matrix<Matrix, SymmGroup>& reduce_column(size_t b2){
            if(!e2[b2]) return *GRID(0,b2);

            std::vector< std::pair< block_matrix<Matrix, SymmGroup>*, int > > rvector;
            for(int b1 = 0; b1 < size1; b1++){
                if(GRID(b1,b2) == NULL) continue;
                rvector.push_back(std::make_pair(GRID(b1,b2), b1));
                GRID(b1,b2) = NULL;
            }
            GRID(0,b2) = ambient::reduce(rvector, [](std::pair<block_matrix<Matrix, SymmGroup>*, int>& dst_pair, 
                                                     std::pair<block_matrix<Matrix, SymmGroup>*, int>& src_pair){
                                                         block_matrix<Matrix, SymmGroup>& dst = *dst_pair.first;
                                                         block_matrix<Matrix, SymmGroup>& src = *src_pair.first;
                                                         ambient::get_scope().set(dst_pair.second);
                                                         for(size_t k = 0; k < src.n_blocks(); ++k)
                                                         dst.match_and_add_block(src[k],
                                                                                 src.left_basis()[k].first, 
                                                                                 src.right_basis()[k].first);
                                                     }).first;
            e2[b2] = false;
            for(int i = 1; i < rvector.size(); i++) delete rvector[i].first;
            return *GRID(0,b2);
        }
        block_matrix<Matrix, SymmGroup>& reduce(){
            std::vector< std::pair< block_matrix<Matrix, SymmGroup>*, int > > rvector;
            std::vector< std::pair< block_matrix<Matrix, SymmGroup>*, int >* > rvector_global;
        
            for(int b2 = 0; b2 < size2; b2++)
            for(int b1 = 0; b1 < size1; b1++){
                if(GRID(b1,b2) == NULL) continue;
                int owner = e2[b2] ? b1 : ambient::scope::permute(b2,mpo.placement_r);
                rvector.push_back(std::make_pair(GRID(b1,b2), owner % ambient::num_workers()));
            }
        
            std::sort(rvector.begin(), rvector.end(), [](const std::pair<block_matrix<Matrix, SymmGroup>*, int>& a, 
                                                         const std::pair<block_matrix<Matrix, SymmGroup>*, int>& b){ 
                                                             return a.second < b.second; 
                                                         });
            int i = 0;
            while(i < rvector.size()){
                int owner = rvector[i].second;
                rvector_global.push_back(&rvector[i++]);
                select_proc(owner);
                while(rvector[i].second == owner && i < rvector.size()){
                    for(size_t k = 0; k < rvector[i].first->n_blocks(); ++k)
                    rvector_global.back()->first->match_and_add_block((*rvector[i].first)[k], 
                                                                      rvector[i].first->left_basis()[k].first, 
                                                                      rvector[i].first->right_basis()[k].first);
                    i++;
                }
            }
            return *ambient::reduce(rvector_global, [](std::pair< block_matrix<Matrix, SymmGroup>*, int >* dst_pair, 
                                                       std::pair< block_matrix<Matrix, SymmGroup>*, int >* src_pair){
                                                           block_matrix<Matrix, SymmGroup>& dst = *dst_pair->first;
                                                           block_matrix<Matrix, SymmGroup>& src = *src_pair->first;
                                                           for(size_t k = 0; k < src.n_blocks(); ++k){
                                                               select_proc(ambient::scope::balance(k,src.n_blocks()));
                                                               dst.match_and_add_block(src[k], 
                                                                                       src.left_basis()[k].first, 
                                                                                       src.right_basis()[k].first);
                                                           }
                                                       })->first;
        }
        MPOTensor<Matrix, SymmGroup> const & mpo;
        mutable std::vector< block_matrix<Matrix, SymmGroup>* > grid;
        std::vector<bool> e2;
        size_t size1;
        size_t size2;
        #undef GRID
    };
}

#endif
