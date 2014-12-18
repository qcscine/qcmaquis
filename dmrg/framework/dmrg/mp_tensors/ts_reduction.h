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

#include <vector>
#include <utility>
#include <alps/numeric/real.hpp>

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/mp_tensors/contractions/non-abelian/gsl_coupling.h"

namespace ts_reduction {
    
    /*
    template<class Matrix, class SymmGroup>
    void reshape_both_to_left(Index<SymmGroup> const & physical_i_left,
                              Index<SymmGroup> const & physical_i_right,
                              Index<SymmGroup> const & left_i,
                              Index<SymmGroup> const & right_i,
                              block_matrix<Matrix, SymmGroup> const & m1,
                              block_matrix<Matrix, SymmGroup> & m2)
    {   
        
        m2 = block_matrix<Matrix, SymmGroup>();
        
        typedef std::size_t size_t;
        typedef typename SymmGroup::charge charge;
        
        Index<SymmGroup> phys2_i = physical_i_left*physical_i_right;
        ProductBasis<SymmGroup> phys_pb(physical_i_left, physical_i_right);
        ProductBasis<SymmGroup> in_left(physical_i_left, left_i);
        ProductBasis<SymmGroup> in_right(physical_i_right, right_i,
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                             -boost::lambda::_1, boost::lambda::_2));
        ProductBasis<SymmGroup> out_left(phys2_i, left_i);
        
       
        for (size_t block = 0; block < m1.n_blocks(); ++block)
        {
            for (size_t s1 = 0; s1 < physical_i_left.size(); ++s1)
            {
                size_t l = left_i.position(SymmGroup::fuse(m1.basis().left_charge(block),
                                                          -physical_i_left[s1].first));
                if(l == left_i.size()) continue;
                for (size_t s2 = 0; s2 < physical_i_right.size(); ++s2)
                {
                    size_t r = right_i.position(SymmGroup::fuse(m1.basis().right_charge(block),
                                                                physical_i_right[s2].first));
                    if(r == right_i.size()) continue;
                    
                    {                            
                        charge s_charge = SymmGroup::fuse(physical_i_left[s1].first, physical_i_right[s2].first);
                        
                        charge out_l_charge = SymmGroup::fuse(s_charge, left_i[l].first);
                        charge out_r_charge = right_i[r].first;
                        
                        //if (! m1.has_block(in_l_charge, in_r_charge) ) continue;
                        
                        // Why is this supposed to work?
                        if (! m2.has_block(out_l_charge, out_r_charge) )
                            m2.insert_block(new Matrix(out_left.size(s_charge, left_i[l].first), right_i[r].second, 0),
                                            out_l_charge, out_r_charge);
                        
                        //maquis::dmrg::detail::reshape_b2l( m2(out_l_charge, out_r_charge), m1(in_l_charge, in_r_charge), 
                        maquis::dmrg::detail::reshape_b2l( m2(out_l_charge, out_r_charge), m1[block], 
                                                           in_left(physical_i_left[s1].first, left_i[l].first), in_right(physical_i_right[s2].first, right_i[r].first),
                                                           out_left(s_charge, left_i[l].first), phys_pb(physical_i_left[s1].first, physical_i_right[s2].first),
                                                           physical_i_left[s1].second, physical_i_right[s2].second, left_i[l].second, right_i[r].second );
                    }
                }
            }
        }
        
    }
    
    template<class Matrix, class SymmGroup>
    void reshape_left_to_both(Index<SymmGroup> const & physical_i_left,
                              Index<SymmGroup> const & physical_i_right,
                              Index<SymmGroup> const & left_i,
                              Index<SymmGroup> const & right_i,
                              block_matrix<Matrix, SymmGroup> const & m1,
                              block_matrix<Matrix, SymmGroup> & m2)
    {   
        
        m2 = block_matrix<Matrix, SymmGroup>();
        
        typedef std::size_t size_t;
        typedef typename SymmGroup::charge charge;
        
        Index<SymmGroup> phys2_i = physical_i_left*physical_i_right;
        ProductBasis<SymmGroup> phys_pb(physical_i_left, physical_i_right);
        ProductBasis<SymmGroup> in_left(phys2_i, left_i);
        
        ProductBasis<SymmGroup> out_right(physical_i_right, right_i,
                                          boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                              -boost::lambda::_1, boost::lambda::_2));
        ProductBasis<SymmGroup> out_left(physical_i_left, left_i);
       
        for (size_t block = 0; block < m1.n_blocks(); ++block)
        {
            size_t r = right_i.position(m1.basis().right_charge(block));
            if(r == right_i.size()) continue;

            for (size_t s1 = 0; s1 < physical_i_left.size(); ++s1)
                for (size_t s2 = 0; s2 < physical_i_right.size(); ++s2)
                {
                    charge s_charge = SymmGroup::fuse(physical_i_left[s1].first, physical_i_right[s2].first);
                    //size_t s_out = phys2_i.position(s_charge);
                    
                    size_t l = left_i.position(SymmGroup::fuse(m1.basis().left_charge(block), -s_charge));
                    if(l == left_i.size()) continue;

                    {
                        charge out_l_charge = SymmGroup::fuse(physical_i_left[s1].first, left_i[l].first);
                        charge out_r_charge = SymmGroup::fuse(-physical_i_right[s2].first, right_i[r].first);
                        //charge in_l_charge = SymmGroup::fuse(s_charge, left_i[l].first);
                        //charge in_r_charge = right_i[r].first;
                        
                        //if (! m1.has_block(in_l_charge, in_r_charge) ) continue;
                        
                        if (! m2.has_block(out_l_charge, out_r_charge) )
                            m2.insert_block(new Matrix(out_left.size(physical_i_left[s1].first, left_i[l].first),
                                                   out_right.size(-physical_i_right[s2].first, right_i[r].first), 0),
                                            out_l_charge, out_r_charge);
                        
                        //maquis::dmrg::detail::reshape_l2b( m2(out_l_charge, out_r_charge), m1(in_l_charge, in_r_charge),
                        maquis::dmrg::detail::reshape_l2b( m2(out_l_charge, out_r_charge), m1[block],
                                                           in_left(s_charge, left_i[l].first), phys_pb(physical_i_left[s1].first, physical_i_right[s2].first),
                                                           out_left(physical_i_left[s1].first, left_i[l].first), out_right(physical_i_right[s2].first, right_i[r].first),
                                                           physical_i_left[s1].second, physical_i_right[s2].second, left_i[l].second, right_i[r].second );
                    }
                }
        }
        
    }
    */

    template <class T, class SymmGroup>
    struct print_values
    {
        void operator()(block_matrix<alps::numeric::matrix<T>, SymmGroup> const & m2) {}
    };

    template <class SymmGroup>
    struct print_values<double, SymmGroup>
    {
        void operator()(block_matrix<alps::numeric::matrix<double>, SymmGroup> const & m2) {
            std::vector<double> vcopy;
            for(int b = 0; b < m2.n_blocks(); ++b)
            {
                std::copy(m2[b].elements().first, m2[b].elements().second, std::back_inserter(vcopy));
            }
            std::sort(vcopy.begin(), vcopy.end());
            maquis::cout << "Number of elements: " << vcopy.size() << std::endl;
            std::copy(vcopy.begin(), vcopy.end(), std::ostream_iterator<double>(maquis::cout, ", "));
            maquis::cout << std::endl;
        }
    };
    
    template<class Matrix, class SymmGroup>
    Index<SymmGroup> reduce_to_right(Index<SymmGroup> const & physical_i_left,
                                     Index<SymmGroup> const & physical_i_right,
                                     Index<SymmGroup> const & left_i,
                                     Index<SymmGroup> const & right_i,
                                     block_matrix<Matrix, SymmGroup> const & m1,
                                     block_matrix<Matrix, SymmGroup> & m2)
    {
        m2 = block_matrix<Matrix, SymmGroup>();
        
        typedef std::size_t size_t;
        typedef typename SymmGroup::charge charge;
        typedef typename Matrix::value_type value_type;
        
        Index<SymmGroup> phys2_i = physical_i_left*physical_i_right;
        ProductBasis<SymmGroup> phys_pb(physical_i_left, physical_i_right);
        ProductBasis<SymmGroup> in_right(phys2_i, right_i,
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                             -boost::lambda::_1, boost::lambda::_2));

        //std::transform(phys_out.begin(), phys_out.end(), phys_out.begin(),
        //               boost::lambda::bind(&std::make_pair<charge, size_t>,
        //               boost::lambda::bind(&std::pair<charge, size_t>::first, boost::lambda::_1), 1) );
        
        maquis::cout << "physical_i_left: " << physical_i_left << std::endl;
        maquis::cout << "physical_i_right: " << physical_i_right << std::endl;
        maquis::cout << "left_i: " << left_i << std::endl;
        maquis::cout << "right_i: " << right_i << std::endl;
        maquis::cout << "phys2_i: " << phys2_i << std::endl;

        for (size_t block = 0; block < m1.n_blocks(); ++block)
        {
            size_t l = left_i.position(m1.basis().left_charge(block));
            assert(l != left_i.size());
            charge lc = left_i[l].first;
            maquis::cout << "lcharge: " << lc << " right_size: " << m1.basis().right_size(block) << std::endl;

            charge in_l_charge = lc;
            charge in_r_charge = m1.basis().right_charge(block);
            
            maquis::cout << "inserting " << left_i[l].second << "x" << m1.basis().right_size(block) << " " << in_l_charge << in_r_charge << std::endl;

            size_t o = m2.insert_block(new Matrix(left_i[l].second, m1.basis().right_size(block), 0) , in_l_charge, in_r_charge);
            Matrix const & in_block = m1[block];
            Matrix & out_block = m2[o];

            for (size_t s = 0; s < phys2_i.size(); ++s)
            {
                charge s_charge = phys2_i[s].first;
                maquis::cout << "  s_charge: " << s_charge << std::endl;;
                size_t r = right_i.position(SymmGroup::fuse(m1.basis().right_charge(block), s_charge));
                if(r == right_i.size()) continue;

                size_t in_right_offset  = in_right(s_charge, right_i[r].first);

                for (size_t s1 = 0; s1 < physical_i_left.size(); ++s1)
                {
                    for (size_t s2 = 0; s2 < physical_i_right.size(); ++s2)
                    {
                        charge phys_c1 = physical_i_left[s1].first, phys_c2 = physical_i_right[s2].first;
                        if (s_charge != SymmGroup::fuse(phys_c1, phys_c2)) continue;
                        
                        size_t in_phys_offset = phys_pb(phys_c1, phys_c2);

                        int jl,jm,jr,S1,S2;
                        S1 = std::abs(phys_c1[1]); S2 = std::abs(phys_c2[1]);
                        jl = in_l_charge[1]; jm = lc[1] + phys_c1[1]; jr = right_i[r].first[1];
                        if (jm < 0) continue;

                        int shift = 0;
                        if ( (jl == jr) && (jl > 0) && (S1 == 1) && (S2 == 1) ) {
                            if (phys_c1[1] < 0)
                               shift = -1; 

                            int j = std::abs(S1-S2);
                            for (int j = std::abs(S1-S2); j <= std::abs(S1+S2); j+=2) {
                                value_type coupling_coeff = std::sqrt(j+1) * std::sqrt(jm+1) * gsl_sf_coupling_6j(jl,jr,j,S2,S1,jm);
                                coupling_coeff = (((jl+jr+S1+S2)/2)%2) ? -coupling_coeff : coupling_coeff;
                                maquis::dmrg::detail::reduce_r(out_block, in_block, coupling_coeff,
                                                               in_right_offset, in_phys_offset, shift + j/2,
                                                               physical_i_left[s1].second, physical_i_right[s2].second,
                                                               left_i[l].second, right_i[r].second);
                            }
                        }
                        else {
                            int j = std::abs(phys_c1[1] + phys_c2[1]);
                            value_type coupling_coeff = std::sqrt(j+1) * std::sqrt(jm+1) * gsl_sf_coupling_6j(jl,jr,j,S2,S1,jm);
                            coupling_coeff = (((jl+jr+S1+S2)/2)%2) ? -coupling_coeff : coupling_coeff;
                            maquis::dmrg::detail::reduce_r(out_block, in_block, coupling_coeff,
                                                           in_right_offset, in_phys_offset, 0,
                                                           physical_i_left[s1].second, physical_i_right[s2].second,
                                                           left_i[l].second, right_i[r].second);
                        }

                        //for (int j = minj; j <= maxj; j+=2)
                        //{
                        //    maquis::cout << "    phys_c1, phys_c2 " << phys_c1 << phys_c2
                        //                 << " access " << in_right_offset << "+" << in_phys_offset*right_i[r].second << "--" << right_i[r].second
                        //                 << "  " << jl<<jr<<j<<S2<<S1<<jm << "  * ";

                        //    value_type coupling_coeff = std::sqrt(j+1) * std::sqrt(jm+1) * gsl_sf_coupling_6j(jl,jr,j,S2,S1,jm);
                        //    coupling_coeff = (((jl+jr+S1+S2)/2)%2) ? -coupling_coeff : coupling_coeff;

                        //    maquis::cout << coupling_coeff << " (phase: " << ((((jl+jr+S1+S2)/2)%2)?-1:1) << ") "
                        //                 << (in_phys_offset + shift + ((use_offset) ? j/2 : 0))*right_i[r].second << "  ";
                        //    //std::copy(in_block.elements().first + in_right_offset + in_phys_offset*right_i[r].second,
                        //    //          in_block.elements().first + in_right_offset + (in_phys_offset+1)*right_i[r].second,
                        //    //          std::ostream_iterator<value_type>(maquis::cout, " "));
                        //    maquis::cout << std::endl;

                        //    maquis::dmrg::detail::reduce_r(out_block, in_block, coupling_coeff,
                        //                                   in_right_offset, in_phys_offset, shift + ((use_offset) ? j/2 : 0),
                        //                                   physical_i_left[s1].second, physical_i_right[s2].second,
                        //                                   left_i[l].second, right_i[r].second);
                        //} // j
                    } // S2
                } // S1
            } // phys2
        } // m1 block

        maquis::cout << std::endl;
        print_values<typename Matrix::value_type, SymmGroup> p;
        //p(m1);
        maquis::cout << std::endl;
        p(m2);

        return phys2_i;
    }
    
    /*
    template<class Matrix, class SymmGroup>
    void reshape_left_to_right(Index<SymmGroup> physical_i,
                               Index<SymmGroup> left_i,
                               Index<SymmGroup> right_i,
                               block_matrix<Matrix, SymmGroup> const & m1,
                               block_matrix<Matrix, SymmGroup> & m2)
    {   }
    */
    
    template<class Matrix, class SymmGroup>
    void reshape_both_to_right(Index<SymmGroup> const & physical_i_left,
                               Index<SymmGroup> const & physical_i_right,
                               Index<SymmGroup> const & left_i,
                               Index<SymmGroup> const & right_i,
                               block_matrix<Matrix, SymmGroup> const & m1,
                               block_matrix<Matrix, SymmGroup> & m2)
    {   
        
        m2 = block_matrix<Matrix, SymmGroup>();
        
        typedef std::size_t size_t;
        typedef typename SymmGroup::charge charge;
        
        Index<SymmGroup> phys2_i = physical_i_left*physical_i_right;
        ProductBasis<SymmGroup> phys_pb(physical_i_left, physical_i_right);
        ProductBasis<SymmGroup> in_left(physical_i_left, left_i);
        ProductBasis<SymmGroup> in_right(physical_i_right, right_i,
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                             -boost::lambda::_1, boost::lambda::_2));
        ProductBasis<SymmGroup> out_right(phys2_i, right_i,
                                          boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                              -boost::lambda::_1, boost::lambda::_2));
       
        for (size_t block = 0; block < m1.n_blocks(); ++block)
        {
            for (size_t s1 = 0; s1 < physical_i_left.size(); ++s1)
            {
                size_t l = left_i.position(SymmGroup::fuse(m1.basis().left_charge(block),
                                                               -physical_i_left[s1].first));
                if(l == left_i.size()) continue;

                for (size_t s2 = 0; s2 < physical_i_right.size(); ++s2)
                {
                    size_t r = right_i.position(SymmGroup::fuse(m1.basis().right_charge(block),
                                                                physical_i_right[s2].first));
                    if(r == right_i.size()) continue;

                    {
                        charge s_charge = SymmGroup::fuse(physical_i_left[s1].first, physical_i_right[s2].first);
                        
                        charge out_l_charge = left_i[l].first;
                        charge out_r_charge = SymmGroup::fuse(-s_charge, right_i[r].first);
                        
                        //if (! m1.has_block(in_l_charge, in_r_charge) ) continue;
                        
                        if (! m2.has_block(out_l_charge, out_r_charge) )
                            m2.insert_block(new Matrix(left_i[l].second, out_right.size(-s_charge, right_i[r].first), 0),
                                            out_l_charge, out_r_charge);
                        
                        maquis::dmrg::detail::reshape_b2r( m2(out_l_charge, out_r_charge), m1[block], 
                                                           in_left(physical_i_left[s1].first, left_i[l].first), in_right(physical_i_right[s2].first, right_i[r].first),
                                                           out_right(s_charge, right_i[r].first), phys_pb(physical_i_left[s1].first, physical_i_right[s2].first),
                                                           physical_i_left[s1].second, physical_i_right[s2].second, left_i[l].second, right_i[r].second );
                    }
                }
            }
        }
        
    }
    
} // namespace ts_reduction
