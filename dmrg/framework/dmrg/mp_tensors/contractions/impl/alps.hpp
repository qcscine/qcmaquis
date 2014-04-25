/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef CONTRACTIONS_IMPL_ALPS_HPP
#define CONTRACTIONS_IMPL_ALPS_HPP

namespace contraction {

    template<class Matrix, class SymmGroup>
    class ContractionGrid {
    public:
        ContractionGrid(MPOTensor<Matrix, SymmGroup> const & mpo, size_t s1, size_t s2)
        : grid(s2), size(s2) {
        }
        void hint(const std::vector<block_matrix<Matrix, SymmGroup> >& t){
        }
        block_matrix<Matrix, SymmGroup>& operator()(size_t b1, size_t b2){
            return grid[b2];
        }
        void multiply_column(size_t b2, const block_matrix<Matrix, SymmGroup>& rhs){
            block_matrix<Matrix, SymmGroup> tmp;
            gemm(grid[b2], rhs, tmp);
            grid[b2] = tmp;
        }
        void multiply_column_trans(size_t b2, const block_matrix<Matrix, SymmGroup>& rhs){
            block_matrix<Matrix, SymmGroup> tmp;
            gemm(transpose(grid[b2]), rhs, tmp);
            grid[b2] = tmp;
        }
        block_matrix<Matrix, SymmGroup> reduce(){
            block_matrix<Matrix, SymmGroup> ret;
            for(int b2 = 0; b2 < size; b2++)
            for(int k = 0; k < grid[b2].n_blocks(); ++k)
                ret.match_and_add_block(grid[b2][k], grid[b2].left_basis()[k].first, grid[b2].right_basis()[k].first);
            return ret;
        }
        Boundary<Matrix, SymmGroup> make_boundary(){
            Boundary<Matrix, SymmGroup> ret;
            ret.resize(size);
            for(int b2 = 0; b2 < size; b2++) ret[b2] = grid[b2];
            return ret;
        }
        mutable std::vector< block_matrix<Matrix, SymmGroup> > grid;
        size_t size;
    };
}

#endif
