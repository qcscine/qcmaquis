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

#ifndef CONTRACTIONS_IMPL_MEMSAVE_HPP
#define CONTRACTIONS_IMPL_MEMSAVE_HPP

namespace contraction {

    template<class Matrix, class SymmGroup>
    class ContractionGrid {
    public:
        ContractionGrid(MPOTensor<Matrix, SymmGroup> const & mpo, size_t s1, size_t s2){
        }
        block_matrix<Matrix, SymmGroup>& operator()(size_t b1, size_t b2){
            return data_;
        }
        void hint(const std::vector<block_matrix<Matrix, SymmGroup> >& t){
            throw std::runtime_error("ContractionGrid::hint not implemented\n");
        }
        void multiply_column(size_t b2, const block_matrix<Matrix, SymmGroup>& rhs){
            throw std::runtime_error("ContractionGrid::multiply_column not implemented\n");
        }
        void multiply_column_trans(size_t b2, const block_matrix<Matrix, SymmGroup>& rhs){
            throw std::runtime_error("ContractionGrid::multiply_column_trans not implemented\n");
        }
        block_matrix<Matrix, SymmGroup> reduce(){
            throw std::runtime_error("ContractionGrid::reduce not implemented\n");
        }
        Boundary<Matrix, SymmGroup> make_boundary(){
            throw std::runtime_error("ContractionGrid::make_boundary not implemented\n");
        }
        mutable block_matrix<Matrix, SymmGroup> data_;
    };
}

#endif
