
#ifndef TENSOR_MATRIX_ITERABLE
#define TENSOR_MATRIX_ITERABLE

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "dmrg/utils/iterator_blas1.h"

namespace detail {

    template<class Matrix, class SymmGroup>
    struct iterable_matrix_impl{
         static void reshape_right_to_left_impl(Index<SymmGroup> & physical_i,
                                                Index<SymmGroup> const & left_i,
                                                Index<SymmGroup> const & right_i,
                                                std::size_t in_right_offset,
                                                std::size_t out_left_offset,
                                                std::size_t s,
                                                std::size_t r,
                                                std::size_t l,
                                                Matrix const & in_block,
                                                Matrix & out_block)
         {
             for (size_t ss = 0; ss < physical_i[s].second; ++ss)
                 for (size_t rr = 0; rr < right_i[r].second; ++rr)
                     // memcpy(&out_block(out_left_offset + ss*left_i[l].second, rr),
                     //        &in_block(0, in_right_offset + ss*right_i[r].second+rr),
                     //        sizeof(typename Matrix::value_type) * left_i[l].second);

//                     Original version,   
                    for(size_t ll = 0; ll < left_i[l].second; ++ll)
                          out_block(out_left_offset + ss*left_i[l].second+ll, rr) = 
                          in_block(ll, in_right_offset + ss*right_i[r].second+rr);
         }

         static void reshape_left_to_right_impl(Index<SymmGroup> & physical_i,
                                                Index<SymmGroup> const & left_i,
                                                Index<SymmGroup> const & right_i,
                                                std::size_t in_left_offset,
                                                std::size_t out_right_offset,
                                                std::size_t s,
                                                std::size_t r,
                                                std::size_t l,
                                                Matrix const & in_block,
                                                Matrix & out_block)
         {
             for (size_t ss = 0; ss < physical_i[s].second; ++ss)
                 for (size_t rr = 0; rr < right_i[r].second; ++rr)
                     for (size_t ll = 0; ll < left_i[l].second; ++ll)
                         out_block(ll, out_right_offset + ss*right_i[r].second+rr) = in_block(in_left_offset + ss*left_i[l].second+ll, rr);
         }


        static void left_boundary_tensor_mpo_impl(Index<SymmGroup> & physical_i,
                                                  Index<SymmGroup> const & left_i,  
                                                  Index<SymmGroup> const & right_i,  
                                                  std::size_t in_left_offset,
                                                  std::size_t out_left_offset,
                                                  std::size_t s1,
                                                  std::size_t s2,
                                                  std::size_t r,
                                                  std::size_t l,
                                                  Matrix const & wblock,
                                                  Matrix const & iblock,
                                                  Matrix & oblock)
        {
            for(size_t ss1 = 0; ss1 < physical_i[s1].second; ++ss1)
                for(size_t ss2 = 0; ss2 < physical_i[s2].second; ++ss2) {
                    typename Matrix::value_type wblock_t = wblock(ss1, ss2);
                    for(size_t rr = 0; rr < right_i[r].second; ++rr) {
                        iterator_axpy(&iblock(in_left_offset + ss1*left_i[l].second, rr),
                                      &iblock(in_left_offset + ss1*left_i[l].second, rr) + left_i[l].second, // bugbug
                                      &oblock(out_left_offset + ss2*left_i[l].second, rr),
                                      wblock_t);
                    }
                }
        }

        static void right_boundary_tensor_mpo_impl(Index<SymmGroup> & physical_i,
                                                   Index<SymmGroup> const & left_i,  
                                                   Index<SymmGroup> const & right_i,  
                                                   std::size_t in_right_offset,
                                                   std::size_t out_right_offset,
                                                   std::size_t s1,
                                                   std::size_t s2,
                                                   std::size_t r,
                                                   std::size_t l,
                                                   Matrix const & wblock,
                                                   Matrix const & iblock,
                                                   Matrix & oblock)
        {
            for (size_t ss1 = 0; ss1 < physical_i[s1].second; ++ss1)
                for (size_t ss2 = 0; ss2 < physical_i[s2].second; ++ss2) {
                    typename Matrix::value_type wblock_t = wblock(ss1, ss2);
                    for (size_t rr = 0; rr < right_i[r].second; ++rr)
                        for (size_t ll = 0; ll < left_i[l].second; ++ll) {
                            oblock(ll, out_right_offset + ss2*right_i[r].second+rr) +=
                            iblock(ll, in_right_offset + ss1*right_i[r].second+rr) * wblock_t;
                        }
                }
        }
 
        static void scalar_norm_impl(Matrix const& M, typename Matrix::value_type & ret)
        {
            using utils::conj;
            for (std::size_t c = 0; c < num_cols(M); ++c)
                for (std::size_t r = 0; r < num_rows(M); ++r)
                    ret += conj(M(r,c)) * M(r,c);
        }

        static void scalar_norm_impl(Matrix const& M1, Matrix const& M2, typename Matrix::value_type & ret)
        {
            using utils::conj;
            for (std::size_t c = 0; c < num_cols(M1); ++c)
                for (std::size_t r = 0; r < num_rows(M1); ++r)
                    ret += conj(M1(r,c)) * M2(r,c);
        }

        static void caculate_bond_renyi_entropies_impl(maquis::types::diagonal_matrix<typename Matrix::value_type> & M, typename maquis::types::associated_real_vector<Matrix>::type& sv) 
        {
            for (typename maquis::types::associated_diagonal_matrix<Matrix>::type::const_element_iterator it = elements(M).first;
                 it != elements(M).second; ++it)
            {
                double a = std::abs(*it);
                if (a > 1e-10)
                    sv.push_back(a*a);
            }
        }
  
        static void left_right_boundary_init_impl(Matrix & M)
        {
//            memset((void*)&M(0,0),1,num_rows(M)*num_cols(M)*sizeof(typename Matrix::value_type));
            for_each(elements(M).first,elements(M).second, boost::lambda::_1 = 1); // boost::lambda ^^' because iterable matrix concept 
        }

    }; // end structure
} // end namespace
#endif
