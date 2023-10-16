/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef CONTRACTIONS_IMPL_ALPS_HPP
#define CONTRACTIONS_IMPL_ALPS_HPP

namespace contraction {

    template<class Matrix, class SymmGroup>
    class ContractionGrid {
    public:
        ContractionGrid(MPOTensor<Matrix, SymmGroup> const & mpo, size_t s1, size_t s2)
        : grid(s2), size(s2), granularity(1) {
        }
        void index_sizes(size_t b2){
            grid[b2].index_sizes();
        }
        static void iterate_reduction_layout(int, int){
        }
        int where(size_t b1, size_t b2){
            return 0;
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
                ret.match_and_add_block(grid[b2][k], grid[b2].basis().left_charge(k), grid[b2].basis().right_charge(k)[k]);
            return ret;
        }
        Boundary<Matrix, SymmGroup> make_boundary(){
            Boundary<Matrix, SymmGroup> ret;
            ret.resize(size);
            for(int b2 = 0; b2 < size; b2++) ret[b2] = grid[b2];
            return ret;
        }
        mutable std::vector< block_matrix<Matrix, SymmGroup> > grid;
        int granularity;
        size_t size;
    };
}

#endif
