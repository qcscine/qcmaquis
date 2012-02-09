
#ifndef P_DENSE_MATRIX_BLOCK_MATRIX_CALLBACKS
#define P_DENSE_MATRIX_BLOCK_MATRIX_CALLBACKS

#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/p_diagonal_matrix.h"
#include "types/p_dense_matrix/algorithms.hpp"

#include "types/utils/matrix_cast.h"

// block_matrix/detail/iterable_matrix.h should be loaded!
// not hard-coding include since it's part of another library
// template<class Matrix, class SymmGroup>
// class block_matrix;


namespace detail {

    // for the block matrix_algorithms
    template <class T, class SymmGroup>
    struct iteretable_diag_impl<maquis::types::p_diagonal_matrix<T>, SymmGroup> {
        static void solver_truncate_impl_zero(block_matrix<maquis::types::p_diagonal_matrix<T>, SymmGroup> const & evals,
                                              size_t Mmax, double cutoff,
                                              size_t* keeps, double & truncated_weight, double & smallest_ev)
        {
            block_matrix<maquis::types::diagonal_matrix<T>, SymmGroup> serial_evals(evals.left_basis(), evals.right_basis());
            for(size_t k=0; k<evals.n_blocks(); ++k){
                serial_evals[k] = maquis::traits::matrix_cast<maquis::types::diagonal_matrix<T> >(evals[k]); // C - a kernel with the vector and the cast inside will be very more efficient ! To DO (only one playout !)
            }
          iteretable_diag_impl<maquis::types::diagonal_matrix<T>, SymmGroup>::solver_truncate_impl_zero(serial_evals, Mmax, cutoff, keeps, truncated_weight, smallest_ev);
        }
    }; // end structure

    //reshape algorithm
    template<class T, class SymmGroup>
    struct iterable_matrix_impl<maquis::types::p_dense_matrix<T>, SymmGroup> {
        // for the reshape right to left
        static void reshape_right_to_left_impl(Index<SymmGroup> & physical_i,
                                               Index<SymmGroup> const & left_i,
                                               Index<SymmGroup> const & right_i,
                                               std::size_t in_right_offset,
                                               std::size_t out_left_offset,
                                               std::size_t s,
                                               std::size_t r,
                                               std::size_t l,
                                               maquis::types::p_dense_matrix<T> const & in_block,
                                               maquis::types::p_dense_matrix<T> & out_block)
        {
            ambient::push(ambient::reshape_r2l_l<T>, ambient::reshape_r2l_c<T>, out_block, in_block, 
                          out_left_offset, in_right_offset, physical_i[s].second, left_i[l].second, right_i[r].second);
            ambient::playout();
        }

        // for the reshape left to right
        static void reshape_left_to_right_impl(Index<SymmGroup> & physical_i,
                                               Index<SymmGroup> const & left_i,
                                               Index<SymmGroup> const & right_i,
                                               std::size_t in_left_offset,
                                               std::size_t out_right_offset,
                                               std::size_t s,
                                               std::size_t r,
                                               std::size_t l,
                                               maquis::types::p_dense_matrix<T> const & in_block,
                                               maquis::types::p_dense_matrix<T> & out_block)
        {
            ambient::push(ambient::reshape_l2r_l<T>, ambient::reshape_l2r_c<T>, in_block, out_block, 
                          in_left_offset, out_right_offset, physical_i[s].second, left_i[l].second, right_i[r].second);
            ambient::playout();
        }

        //for the contraction of the left boundary tensor
        static void left_boundary_tensor_mpo_impl(Index<SymmGroup> & physical_i,
                                                  Index<SymmGroup> const & left_i,  
                                                  Index<SymmGroup> const & right_i,  
                                                  std::size_t in_left_offset,
                                                  std::size_t out_left_offset,
                                                  std::size_t s1,
                                                  std::size_t s2,
                                                  std::size_t r,
                                                  std::size_t l,
                                                  maquis::types::p_dense_matrix<T> const & wblock,
                                                  maquis::types::p_dense_matrix<T> const & iblock,
                                                  maquis::types::p_dense_matrix<T> & oblock)
        {
            ambient::push(ambient::lb_tensor_mpo_l<T>, ambient::lb_tensor_mpo_c<T>,
                          oblock, iblock, wblock, out_left_offset, in_left_offset, 
                          physical_i[s1].second, physical_i[s2].second, left_i[l].second, right_i[r].second);
            ambient::playout(); 
        }

        //for the contraction of the left boundary tensor
        static void right_boundary_tensor_mpo_impl(Index<SymmGroup> & physical_i,
                                                   Index<SymmGroup> const & left_i,  
                                                   Index<SymmGroup> const & right_i,  
                                                   std::size_t in_right_offset,
                                                   std::size_t out_right_offset,
                                                   std::size_t s1,
                                                   std::size_t s2,
                                                   std::size_t r,
                                                   std::size_t l,
                                                   maquis::types::p_dense_matrix<T> const & wblock,
                                                   maquis::types::p_dense_matrix<T> const & iblock,
                                                   maquis::types::p_dense_matrix<T> & oblock)
        {
            ambient::push(ambient::rb_tensor_mpo_l<T>, ambient::rb_tensor_mpo_c<T>,
                         oblock, iblock, wblock, out_right_offset, in_right_offset, 
                         physical_i[s1].second, physical_i[s2].second, left_i[l].second, right_i[r].second);
            ambient::playout();
        }

        static void scalar_norm_impl(Matrix& M, typename Matrix::value_type & ret){ // not const due to nullcut 
            M.nullcut(); // not counting redunant elements of workgroup
            typename Matrix::value_type* norm = &ret;
            ambient::push(ambient::scalar_norm_l<T>, ambient::scalar_norm_c<T>, M, norm);
            ambient::playout(); // execution weight: 452
        }

        static void scalar_norm_impl(maquis::types::p_dense_matrix<T> & M1, maquis::types::p_dense_matrix<T> & M2, typename Matrix::value_type & ret){ // not const due to nullcut
            M1.nullcut(); // not counting redunant elements of workgroup
            typename Matrix::value_type* overlap = &ret;
            ambient::push(ambient::scalar_overlap_l<T>, ambient::scalar_overlap_c<T>, M1, M2, overlap);
            ambient::playout(); // execution weight: 452
        }

        static void caculate_bond_renyi_entropies_impl(maquis::types::p_diagonal_matrix<typename Matrix::value_type> const& M, std::vector<typename Matrix::value_type>& sv){
//          maquis::types::algorithms::copy_sqr_gt<Matrix>(sv, M, 1e-10);
            std::vector<typename Matrix::value_type>* sc_ptr = &sv;
            double prec(1e-10); 
            ambient::push(ambient::push_back_sqr_gt_l, ambient::push_back_sqr_gt_c, sc_ptr, M.get_data(), prec);
            ambient::playout();
        }

        static void left_right_boundary_init_impl(maquis::types::p_dense_matrix<T> & M){
            double one(1.0);
            ambient::push(ambient::initv_l<T>,ambient::initv_c<T>, M, one);
            ambient::playout();
        }
    }; // end structure
}
#endif
