/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include <vector>
#include <utility>

#include "dmrg/block_matrix/symmetry/gsl_coupling.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

namespace ts_reduction {

    namespace detail {

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
    }

    template<class Matrix, class SymmGroup>
    Index<SymmGroup> reduce_right(Index<SymmGroup> const & physical_i_left,
                                  Index<SymmGroup> const & physical_i_right,
                                  Index<SymmGroup> const & left_i,
                                  Index<SymmGroup> const & right_i,
                                  block_matrix<Matrix, SymmGroup> const & m1,
                                  block_matrix<Matrix, SymmGroup> & m2)
    {
        m2 = block_matrix<Matrix, SymmGroup>();

        typedef std::size_t size_t;
        typedef typename SymmGroup::subcharge spin_t;
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

        for (size_t block = 0; block < m1.n_blocks(); ++block)
        {
            charge lc = m1.basis().left_charge(block);
            size_t left_size = m1.basis().left_size(block);

            charge in_r_charge = m1.basis().right_charge(block);

            size_t o = m2.insert_block(new Matrix(left_size, m1.basis().right_size(block), 0) , lc, in_r_charge);
            Matrix const & in_block = m1[block];
            Matrix & out_block = m2[o];

            for (size_t s = 0; s < phys2_i.size(); ++s)
            {
                charge s_charge = phys2_i[s].first;
                size_t r = right_i.position(SymmGroup::fuse(in_r_charge, s_charge));
                if(r == right_i.size()) continue;

                size_t in_right_offset  = in_right(s_charge, right_i[r].first);
                size_t right_size = right_i[r].second;

                for (size_t s1 = 0; s1 < physical_i_left.size(); ++s1)
                {
                    for (size_t s2 = 0; s2 < physical_i_right.size(); ++s2)
                    {
                        charge phys_c1 = physical_i_left[s1].first, phys_c2 = physical_i_right[s2].first;
                        if (s_charge != SymmGroup::fuse(phys_c1, phys_c2)) continue;

                        size_t in_phys_offset = phys_pb(phys_c1, phys_c2);

                        spin_t jl,jm,jr,S1,S2;
                        S1 = std::abs(SymmGroup::spin(phys_c1));
                        S2 = std::abs(SymmGroup::spin(phys_c2));
                        jl = SymmGroup::spin(lc);
                        jm = SymmGroup::spin(lc) + SymmGroup::spin(phys_c1);
                        jr = SymmGroup::spin(right_i[r].first);

                        if (jm < 0) continue;

                        if ( (jl == jr) && (jl > 0) && (S1 == 1) && (S2 == 1) ) {
                            size_t base_offset = (SymmGroup::spin(phys_c1) == 1) ? in_phys_offset : in_phys_offset - 1;

                            //  MPS right_paired layout                 |    <2,0> MPS sector: (c1×c2)    |             this is a dense matrix
                            //                                          right_size, r_charge = lc + phys_charge              ↓
                            // lc (jl)   -> |---------------------------    ^    ----^--------^-------^------------------------------------|
                            // left_size -> | (onsite charges > <2,0>)  | 20x00 | 11x1-1 | 1-1x11 | 00x20 |   (onsite charges < <2,0>)     |
                            //              |---------------------------------------A---------B--------------------------------------------|
                            //
                            //               in_right_offset----------->
                            //                                         > in_phys_offsets (values 0 to 3 max), multiply by left_size
                            //                                         --------> base_offset (for recoupling A and B)
                            //                                         ------------------>
                            //                                         -------------------------->
                            //
                            //                                  A_reduced(spin = 0) = 6j(..) * A(jm = jl+1) + 6j(..) * B(jm = jl-1)
                            //                                  B_reduced(spin = 1) = 6j(..) * A(jm = jl+1) + 6j(..) * B(jm = jl-1)

                            for (spin_t j = std::abs(S1-S2); j <= std::abs(S1+S2); j+=2) {
                                size_t out_phys_offset = base_offset + j/2;
                                value_type coupling_coeff = std::sqrt((j+1) * (jm+1)) * WignerWrapper::gsl_sf_coupling_6j(jl,jr,j,S2,S1,jm);
                                coupling_coeff = (((jl+jr+S1+S2)/2)%2) ? -coupling_coeff : coupling_coeff;
                                maquis::dmrg::detail::reduce_r(out_block, in_block, coupling_coeff,
                                                               in_right_offset, in_phys_offset, out_phys_offset,
                                                               physical_i_left[s1].second, physical_i_right[s2].second,
                                                               left_size, right_size);
                            }
                        }
                        else {
                            spin_t j = std::abs(SymmGroup::spin(phys_c1) + SymmGroup::spin(phys_c2));
                            value_type coupling_coeff = std::sqrt((j+1) * (jm+1)) * WignerWrapper::gsl_sf_coupling_6j(jl,jr,j,S2,S1,jm);
                            coupling_coeff = (((jl+jr+S1+S2)/2)%2) ? -coupling_coeff : coupling_coeff;
                            maquis::dmrg::detail::reduce_r(out_block, in_block, coupling_coeff,
                                                           in_right_offset, in_phys_offset, in_phys_offset,
                                                           physical_i_left[s1].second, physical_i_right[s2].second,
                                                           left_size, right_size);
                        }

                    } // S2
                } // S1
            } // phys2
        } // m1 block

        return phys2_i;
    }

    template<class Matrix, class SymmGroup>
    Index<SymmGroup> unreduce_left(Index<SymmGroup> const & physical_i_left,
                                   Index<SymmGroup> const & physical_i_right,
                                   Index<SymmGroup> const & left_i,
                                   Index<SymmGroup> const & right_i,
                                   block_matrix<Matrix, SymmGroup> const & m1,
                                   block_matrix<Matrix, SymmGroup> & m2)
    {
        m2 = block_matrix<Matrix, SymmGroup>();

        typedef std::size_t size_t;
        typedef typename SymmGroup::subcharge spin_t;
        typedef typename SymmGroup::charge charge;
        typedef typename Matrix::value_type value_type;

        Index<SymmGroup> phys2_i = physical_i_left*physical_i_right;
        ProductBasis<SymmGroup> phys_pb(physical_i_left, physical_i_right);
        ProductBasis<SymmGroup> in_left(phys2_i, left_i);

        for (size_t block = 0; block < m1.n_blocks(); ++block)
        {
            charge rc = m1.basis().right_charge(block);
            size_t right_size = m1.basis().right_size(block);

            charge in_l_charge = m1.basis().left_charge(block);

            size_t o = m2.insert_block(new Matrix(m1.basis().left_size(block), right_size, 0) , in_l_charge, rc);
            Matrix const & in_block = m1[block];
            Matrix & out_block = m2[o];

            for (size_t s = 0; s < phys2_i.size(); ++s)
            {
                charge s_charge = phys2_i[s].first;
                size_t l = left_i.position(SymmGroup::fuse(in_l_charge, -s_charge));
                if(l == left_i.size()) continue;

                size_t in_left_offset = in_left(s_charge, left_i[l].first);

                for (size_t s1 = 0; s1 < physical_i_left.size(); ++s1)
                {
                    for (size_t s2 = 0; s2 < physical_i_right.size(); ++s2)
                    {
                        charge phys_c1 = physical_i_left[s1].first, phys_c2 = physical_i_right[s2].first;
                        if (s_charge != SymmGroup::fuse(phys_c1, phys_c2)) continue;

                        size_t in_phys_offset = phys_pb(phys_c1, phys_c2);

                        spin_t jl,j,jr,S1,S2;
                        S1 = std::abs(SymmGroup::spin(phys_c1));
                        S2 = std::abs(SymmGroup::spin(phys_c2));
                        jl = SymmGroup::spin(left_i[l].first);
                        jr = SymmGroup::spin(rc);

                        if ( (jl == jr) && (jl > 0) && (S1 == 1) && (S2 == 1) ) {
                            spin_t j = (SymmGroup::spin(phys_c1) == 1) ? 0 : 2;
                            size_t base_offset = (j==0) ? in_phys_offset : in_phys_offset - 1;

                            //  MPS left_paired layout (rotate -90°)     |    <2,0> MPS sector: (c1×c2)    |           this is a dense matrix
                            //                                           left_size                                        ↓
                            //               |---------------------------    ^    ----^--------^-------^------------------------------------|
                            // right_size -> | (onsite charges > <2,0>)  | 20x00 | 11x1-1 | 1-1x11 | 00x20 |   (onsite charges < <2,0>)     |
                            //               |-------------------------------------A(S=0)---B(S=1)------------------------------------------|
                            //
                            //                 in_left_offset----------->
                            //                                           > in_phys_offsets (values 0 to 3 max), multiply by left_size
                            //                                           --------> base_offset (for recoupling A and B)
                            //                                           ------------------>
                            //                                           -------------------------->
                            //
                            //                                  A_unreduced(jm = jl+1) = 6j(..) * A(S=0) + 6j(..) * B(S=2)
                            //                                  B_unreduced(jm = jl-1) = 6j(..) * A(S=0) + 6j(..) * B(S=2)

                            for (spin_t jm = jl - 1; jm <= jl + 1; jm+=2) {
                                size_t out_phys_offset = (jm == jl - 1) ? base_offset + 1 : base_offset;
                                value_type coupling_coeff = std::sqrt((j+1) * (jm+1)) * WignerWrapper::gsl_sf_coupling_6j(jl,jr,j,S2,S1,jm);
                                coupling_coeff = (((jl+jr+S1+S2)/2)%2) ? -coupling_coeff : coupling_coeff;
                                maquis::dmrg::detail::reduce_l(out_block, in_block, coupling_coeff,
                                                               in_left_offset, in_phys_offset, out_phys_offset,
                                                               physical_i_left[s1].second, physical_i_right[s2].second,
                                                               left_i[l].second, right_size);
                            }
                        }
                        else {
                            spin_t j = std::abs(SymmGroup::spin(phys_c1) + SymmGroup::spin(phys_c2));
                            spin_t jm = jl + SymmGroup::spin(phys_c1);
                            if (jm < 0) continue;
                            value_type coupling_coeff = std::sqrt((j+1) * (jm+1)) * WignerWrapper::gsl_sf_coupling_6j(jl,jr,j,S2,S1,jm);
                            coupling_coeff = (((jl+jr+S1+S2)/2)%2) ? -coupling_coeff : coupling_coeff;
                            maquis::dmrg::detail::reduce_l(out_block, in_block, coupling_coeff,
                                                           in_left_offset, in_phys_offset, in_phys_offset,
                                                           physical_i_left[s1].second, physical_i_right[s2].second,
                                                           left_i[l].second, right_size);
                        }

                    } // S2
                } // S1
            } // phys2
        } // m1 block

        return phys2_i;
    }

} // namespace ts_reduction
