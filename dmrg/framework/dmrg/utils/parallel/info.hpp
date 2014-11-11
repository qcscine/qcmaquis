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

#ifndef UTILS_PARALLEL_INFO_HPP
#define UTILS_PARALLEL_INFO_HPP

namespace parallel {

    template<class Matrix, class SymmGroup>
    void print_distribution(const Boundary<Matrix,SymmGroup>& boundary) const {
        if(!parallel::master()) return;
        double total = 0;
        int loop_max = boundary.aux_dim();
        for(int b = 0; b < loop_max; ++b) total += boundary[b].num_elements();
        printf("%.2f GB:", total*sizeof(value_type)/1024/1024/1024);
        for(int p = 0; p < traits::size(); ++p){
            double part = 0;
            for(int b = 0; b < loop_max; ++b)
            for(int i = 0; i < boundary[b].n_blocks(); ++i){
                if(!ambient::weak(boundary[b][i][0]) && traits::placement(boundary[b][i][0]) == *traits::to_iterator(p))
                    part += num_rows(boundary[b][i])*num_cols(boundary[b][i]);
            }
            printf(" %.1f%%", 100*part/total);
        }
        printf("\n");
    }

    template<class Matrix, class SymmGroup>
    void print_distribution(const block_matrix<Matrix, SymmGroup>& m) const {
        if(!parallel::master()) return;
        double total = num_elements();
        printf("%.2f GB:", total*sizeof(typename Matrix::value_type)/1024/1024/1024);
        for(int p = 0; p < traits::size(); ++p){
            double part = 0;
            for(int i = 0; i < m.n_blocks(); ++i){
                if(!ambient::weak(m[i][0]) && traits::placement(m[i][0]) == *traits::to_iterator(p))
                    part += num_rows(m[i])*num_cols(m[i]);
            }
            printf(" %.1f%%", 100*part/total);
        }
        printf("\n");
    }

}

#endif
