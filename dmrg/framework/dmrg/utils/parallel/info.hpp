/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
