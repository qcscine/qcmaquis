/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef REF_H_DIAG_H
#define REF_H_DIAG_H

    //template<class Matrix, class SymmGroup>
    //void ref_diag(SiteProblem<Matrix, SymmGroup> const & H,
    //              MPSTensor<Matrix, SymmGroup> x)
    //{
    //}

    template<class Matrix, class SymmGroup>
    void ref_diag(SiteProblem<Matrix, SymmGroup> const & H,
                  MPSTensor<Matrix, SymmGroup> x)
    {
        typedef typename Matrix::value_type value_type;
        typedef typename SymmGroup::charge charge;

        x.make_left_paired();
        x.multiply_by_scalar(.0);
        block_matrix<Matrix, SymmGroup> & bm = x.data();


        Index<SymmGroup> const & physical_i = x.site_dim(),
                               & left_i = x.row_dim();
        Index<SymmGroup> right_i = x.col_dim(),
                         out_left_i = physical_i * left_i;

        common_subset(out_left_i, right_i);
        ProductBasis<SymmGroup> out_left_pb(physical_i, left_i);
        ProductBasis<SymmGroup> in_right_pb(physical_i, right_i,
                                boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                        -boost::lambda::_1, boost::lambda::_2));

        ///////////////////////////////////////////////////////////////////////
        block_matrix<Matrix, SymmGroup> ret2;
        for (size_t b2 = 0; b2 < H.right.aux_dim(); ++b2)
        {
            block_matrix<Matrix, SymmGroup> la = contraction::SU2::h_diag2(b2, H.left, H.mpo, bm.basis(),
                                                                           x.row_dim(), x.col_dim(), physical_i,
                                                                           in_right_pb, out_left_pb);

            for (size_t block = 0; block < la.n_blocks(); ++block)
            {
                charge in_r_charge = la.basis().right_charge(block);
                size_t rblock = H.right[b2].find_block(in_r_charge, in_r_charge);
                if (rblock != H.right[b2].n_blocks())
                {
                    for (size_t c = 0; c < num_cols(la[block]); ++c)
                        std::transform(la[block].col(c).first, la[block].col(c).second, la[block].col(c).first,
                                       boost::lambda::_1*= H.right[b2][rblock](c,c));
                    ret2.match_and_add_block(la[block], in_r_charge, in_r_charge);
                }
            }
        }
        ///////////////////////////////////////////////////////////////////////

        for (size_t block = 0; block < bm.n_blocks(); ++block)
        {
            size_t r = right_i.position(bm.basis().right_charge(block));
            assert(r == block);
            if(r == right_i.size()) continue;
            charge in_r_charge = right_i[r].first;
            charge in_l_charge = bm.basis().left_charge(block);

            maquis::cout << in_l_charge << in_r_charge << std::endl;

            for (size_t s = 0; s < physical_i.size(); ++s)
            {
                charge phys_charge = physical_i[s].first;
                size_t l = left_i.position(SymmGroup::fuse(in_l_charge, -phys_charge));
                if(l == left_i.size()) continue;

                size_t in_left_offset = out_left_pb(phys_charge, left_i[l].first);

                for (size_t in_phys_offset = 0; in_phys_offset < physical_i[s].second; ++in_phys_offset)
                {
                    for (size_t i = 0; i < left_i[l].second; ++i)
                    for (size_t j = 0; j < num_cols(bm[block]); ++j)
                    {
                        value_type ret = 0.0;
                        for (size_t b2 = 0; b2 < H.right.aux_dim(); ++b2)
                        {
                            value_type la = contraction::SU2::h_diag(b2, H.left, H.mpo,
                                                                     x.row_dim(), x.col_dim(), physical_i,
                                                                     in_right_pb, out_left_pb,
                                                                     left_i[l].first,
                                                                     phys_charge,
                                                                     in_phys_offset,
                                                                     i);
                            size_t rblock = H.right[b2].find_block(in_r_charge, in_r_charge);
                            if (rblock != H.right[b2].n_blocks())
                                ret += la * H.right[b2][rblock](j,j);

                        }

                        size_t total_index = i + in_left_offset + in_phys_offset * left_i[l].second;
                        bm[block](total_index, j) = 1;    
                        MPSTensor<Matrix, SymmGroup> prod;
                        ietl::mult(H, x, prod);
                        value_type ret_ref = prod.data()[block](total_index, j);
                        bm[block](total_index,j) = 0;    

                        size_t ret2bl = ret2.find_block(in_r_charge, in_l_charge);
                        maquis::cout << "  " << total_index << "," << j << "  " << ret_ref << " " << ret2[ret2bl](total_index, j) << std::endl;

                    } // i,j
                } // phys offset
            } // physical_i s
        }

        //for (size_t b = 0; b < bm.n_blocks(); ++b)
        //{
        //    maquis::cout << bm.basis().left_charge(b) << bm.basis().right_charge(b) << std::endl;
        //    for (size_t i = 0; i < num_rows(bm[b]); ++i)
        //    for (size_t j = 0; j < num_cols(bm[b]); ++j)
        //    {
        //        bm[b](i,j) = 1;    
        //        MPSTensor<Matrix, SymmGroup> prod;
        //        ietl::mult(H, x, prod);
        //        maquis::cout << "  " << i << "," << j << "  " << prod.data()[b](i,j) << std::endl;
        //        bm[b](i,j) = 0;    
        //    }
        //}
    }

    //template<class Matrix, class SymmGroup>
    //void exp_diag(SiteProblem<Matrix, SymmGroup> const & H,
    //              MPSTensor<Matrix, SymmGroup> x)
    //{
    //    typedef typename Matrix::value_type value_type;

    //    value_type la = contraction::SU2::h_diag(0, H.left, H.right, 
    //}

#endif 
