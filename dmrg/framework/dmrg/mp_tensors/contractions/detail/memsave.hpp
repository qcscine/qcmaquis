/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef CONTRACTIONS_IMPL_MEMSAVE_HPP
#define CONTRACTIONS_IMPL_MEMSAVE_HPP

namespace contraction {

    template<class Matrix, class SymmGroup>
    class ContractionGrid {
    public:
        ContractionGrid(MPOTensor<Matrix, SymmGroup> const & mpo, size_t s1, size_t s2) : granularity(1) {
        }
        block_matrix<Matrix, SymmGroup>& operator()(size_t b1, size_t b2){
            return data_;
        }
        void index_sizes(size_t){
            data_.index_sizes();
        }
        static void iterate_reduction_layout(int, int){
        }
        int where(size_t b1, size_t b2){
            return 0;
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
        int granularity;
    };
}

#endif
