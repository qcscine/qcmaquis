/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2018         Leon Freitag <lefreita@ethz.ch>
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
#ifndef NRDMDERIVATIVE_H
#define NRDMDERIVATIVE_H

#include "dmrg/models/measurements/tagged_nrankrdm.h"
#include "dmrg/mp_tensors/twositetensor.h"

    namespace measurements
    {
        template <class Matrix, class SymmGroup>
        class NRDMDerivative : public TaggedNRankRDM<Matrix, SymmGroup>
        {
        typedef TaggedNRankRDM<Matrix, SymmGroup> base;
        typedef typename base::tag_vec tag_vec;
        typedef typename base::scaled_bond_term scaled_bond_term;
        typedef typename base::positions_type positions_type;


        using base::ext_labels;
        public:
            // Specialization to call the correct constructor of the base class
            // TODO: Check if this compiles with disabled SU2U1/SU2U1PG!
            // SU2U1: symm_traits::HasSU2<SU2U1PG> is equivalent to boost::true_type
            NRDMDerivative(symm_traits::HasSU2<SU2U1PG>, std::string const& name_, const Lattice & lat,
                       boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_,
                       typename TermMakerSU2<Matrix, SymmGroup>::OperatorCollection const & op_collection_,
                       positions_type const& positions_ = positions_type())
                       : base(name_, lat, tag_handler_, op_collection_, positions_) {};

            // 2U1 and other symmetry groups
            NRDMDerivative(boost::false_type, std::string const& name_, const Lattice & lat,
                       boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_,
                       tag_vec const & identities_, tag_vec const & fillings_, std::vector<scaled_bond_term> const& ops_,
                       bool half_only_, positions_type const& positions_ = positions_type())
                        : base(name_, lat, tag_handler_, identities_, fillings_, ops_, half_only_, positions_) {};




            virtual void evaluate(MPS<Matrix, SymmGroup> const& ket_mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
            {
                if (this->name() == "onerdmderiv")
                    measure_1rdm_derivative(ket_mps);
            }
        protected:
            measurement<Matrix, SymmGroup>* do_clone() const
            {
                return new NRDMDerivative(*this);
            }

            // Calculates RDM derivatives with respect to MPS parameters
            // If twosite is not set, MPS parameters are calculated only from MPSTensor at site [site]
            // otherwise they are calculated from a TwoSiteTensor at [site] and [site+1]

            // TODO: 0 is only temporarily the default argument for site until I figure out how to pass the site reasonably
            void measure_1rdm_derivative(MPS<Matrix, SymmGroup> mps, int site = 0, bool twosite = true)
            {

                MPS<Matrix, SymmGroup> mps_aux = mps;

                if (twosite) // MPS parameters from two-site tensors
                {
                    // forbid the last site because we can't build a two-site tensor
                    if (site > mps.length() - 1)
                        throw std::runtime_error("site > L-1 not allowed in building RDM derivatives with two-site tensors");

                    // Prepare the two-site tensor from two sites of the MPS
                    TwoSiteTensor<Matrix, SymmGroup> tst(mps_aux[site], mps_aux[site+1]);
                    MPSTensor<Matrix, SymmGroup> twin_mps = tst.make_mps();

                    // Zero all elements
                    twin_mps.multiply_by_scalar(0.);

                    // Loop over all elements. Use three indices in the same way Yingjin does:
                    // i runs over the number of blocks
                    // j and k over the rows and column of each block

                    for (int i = 0; i < twin_mps.data().n_blocks(); i++)
                    for (int j = 0; j < twin_mps.data()[i].num_rows(); j++)
                    for (int k = 0; k < twin_mps.data()[i].num_cols(); k++)
                    {
                        // write indices into ext_indices so that they get picked up during measurements
                        this->ext_labels = { i, j, k };

                        // Prepare the auxiliary state: set the corresponding element to 1 in the TwoSiteTensor
                        twin_mps.data()[i](j,k) = 1.0;

                        // Incorporate the modified TwoSiteTensor back into the MPS to form auxiliary state
                        // Note that we provide a m value as large as possible (maxint) and a negative cutoff in order to prevent truncation.

                        truncation_results trunc;
                        boost::tie(mps_aux[site], mps_aux[site+1], trunc) = tst.split_mps_l2r(std::numeric_limits<int>::max(), -1.0);
                        mps_aux[site].make_left_paired();
                        mps_aux[site+1].make_left_paired();
                        // TODO: Check if the resulting aux_mps is correct!

                        // measure the transition 1-RDM <mps_aux|c+c|mps>
                        // TODO: Note that we will need both <mps_aux|c+c|mps> and <mps|c+c|mps_aux> for symmetrised derivatives.
                        // The symmetrisation will be taken care for later
                        this->measure_correlation(mps_aux, mps);

                        // reset the current element back to 0
                        twin_mps.data()[i](j,k) = 0.0;
                    }
                }
            }
        };
    } // measurements

#endif