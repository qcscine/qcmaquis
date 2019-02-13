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
            // symm_traits::HasSU2<SU2U1> and symm_traits::HasSU2<SU2U1PG> yields boost::true_type
            NRDMDerivative(int lr_site_, boost::true_type, std::string const& name_, const Lattice & lat,
                       boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_,
                       typename TermMakerSU2<Matrix, SymmGroup>::OperatorCollection const & op_collection_,
                       positions_type const& positions_ = positions_type())
                       : base(name_, lat, tag_handler_, op_collection_, positions_), lr_site(lr_site_) {};

            // 2U1 and other symmetry groups
            NRDMDerivative(int lr_site_, boost::false_type, std::string const& name_, const Lattice & lat,
                       boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_,
                       tag_vec const & identities_, tag_vec const & fillings_, std::vector<scaled_bond_term> const& ops_,
                       bool half_only_, positions_type const& positions_ = positions_type())
                       : base(name_, lat, tag_handler_, identities_, fillings_, ops_, half_only_, positions_), lr_site(lr_site_) {};

            virtual void evaluate(MPS<Matrix, SymmGroup> const& ket_mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
            {
                // Wrapper for the TDM evaluation function
                std::function<void(MPS<Matrix, SymmGroup> const &, MPS<Matrix, SymmGroup> const &)> RDMEvaluator;

                // Assign the appropriate TDM evaluation function to RDMEvaluator based on which derivative should be calculated
                if (this->name() == "onerdmderivR")
                    RDMEvaluator = std::bind(&NRDMDerivative<Matrix, SymmGroup>::measure_correlation, this, std::placeholders::_1, std::placeholders::_2);

                if (this->name() == "twordmderivR")
                    RDMEvaluator = std::bind(&NRDMDerivative<Matrix, SymmGroup>::measure_2rdm, this, std::placeholders::_1, std::placeholders::_2);

                // for left RDM derivatives, the bra and ket in the RDM evaluation are swapped
                if (this->name() == "onerdmderivL")
                    RDMEvaluator = std::bind(&NRDMDerivative<Matrix, SymmGroup>::measure_correlation, this, std::placeholders::_2, std::placeholders::_1);

                if (this->name() == "twordmderivL")
                    RDMEvaluator = std::bind(&NRDMDerivative<Matrix, SymmGroup>::measure_2rdm, this, std::placeholders::_2, std::placeholders::_1);

                measure_derivative(ket_mps, RDMEvaluator, lr_site);
            }
        protected:
            measurement<Matrix, SymmGroup>* do_clone() const
            {
                return new NRDMDerivative(*this);
            }

            // MEASURE_DERIVATIVE
            //
            // Generic function to calculate RDM derivatives with respect to MPS parameters
            // Constructs an auxiliary state from MPS [mps] (as described in doi:10.1021/acs.jctc.6b01118)
            // and calculates RDM derivatives as TDMs between [mps] and the auxiliary state
            // If [twosite] is not set, MPS parameters are calculated only from an MPSTensor at site [site],
            // otherwise they are calculated from a TwoSiteTensor at [site] and [site+1]
            //
            // The function to calculate RDMs ( for which the derivative should be calculated)
            // is passed as [RDMevaluator], which allows calculating derivatives for all available RDMs

            // TODO: 0 is only temporarily the default argument for site until I figure out how to pass the site reasonably
            void measure_derivative(MPS<Matrix, SymmGroup> const& mps, std::function<void(MPS<Matrix, SymmGroup> const &, MPS<Matrix, SymmGroup> const &)> RDMEvaluator, int site = 0, bool twosite = true)
            {

                MPS<Matrix, SymmGroup> mps_aux = mps;
                // The MPS needs to be canonized up to the site before the measurement to have the same local basis in all the states
                // However, we canonize all MPS simultaneously, so we don't need to canonize here again.
                // mps_aux.canonize(site);

                if (twosite) // MPS parameters from two-site tensors
                {
                    // forbid the last site because we can't build a two-site tensor
                    if (site > mps.length() - 1)
                        throw std::runtime_error("site > L-1 not allowed in building RDM derivatives with two-site tensors");

                    // Prepare the two-site tensor from two sites of the MPS
                    TwoSiteTensor<Matrix, SymmGroup> tst(mps_aux[site], mps_aux[site+1]);

                    // we need to keep the original elements for the sign correction
                    TwoSiteTensor<Matrix, SymmGroup> tst_orig = tst;

                    // To keep the consistency with other measurements (local Hamiltonian), or the yingjin-devel branch
                    // which apparently works with left-paired two site tensors, we introduce left pairing
                    // Note that for some reason, pairing may introduce an additional block with a zero element!
                    tst.make_left_paired();
                    tst_orig.make_left_paired();

                    // Zero all elements
                    tst.data()*= 0.0;

                    // Loop over all elements. Use three indices in the same way Yingjin does:
                    // i runs over the number of blocks
                    // j and k over the rows and column of each block

                    for (int i = 0; i < tst.data().n_blocks(); i++)
                    for (int j = 0; j < tst.data()[i].num_rows(); j++)
                    for (int k = 0; k < tst.data()[i].num_cols(); k++)
                    {
                        // write indices into ext_indices so that they get picked up during measurements
                        this->ext_labels = { i, j, k };

                        // TODO: Check for zero elements in the original MPS to avoid unnecessary calculations

                        // Prepare the auxiliary state: set the corresponding element to 1 in the TwoSiteTensor
                        // Note re sign correction: to maintain the correct sign of the RDM derivatives, the sign of the original TwoSiteTensor element must be retained
                        // This is much more efficient than calculating the overlap between the MPS and the auxiliary MPS and correcting the derivative signs with the sign of the overlap
                        tst.data()[i](j,k) = std::copysign(1.0, maquis::real(tst_orig.data()[i](j,k)));

                        // Incorporate the modified TwoSiteTensor back into the MPS to form auxiliary state
                        // Note that we provide a m value as large as possible (maxint) and a negative cutoff in order to prevent truncation.

                        truncation_results trunc;
                        boost::tie(mps_aux[site], mps_aux[site+1], trunc) = tst.split_mps_l2r(std::numeric_limits<int>::max(), -1.0);

                        // Fix the left pairing since manipulation of the MPS messes it up
                        tst.make_left_paired();

                        // FIXME: !!!there is a horrible race condition wrt pairing somewhere in expval() that needs to be fixed!!!
                        // Currently, the code doesn't run for OMP_NUM_THREADS > 1 !
                        // The lines below alleviate the problem, but not always!
                        // mps_aux[site].make_left_paired();
                        // mps_aux[site+1].make_left_paired();
                        // mps[site].make_left_paired();
                        // mps[site+1].make_left_paired();

                        // measure the transition RDM <mps_aux|c+...c...|mps>
                        // TODO: Note that we will need both <mps_aux|c+...c...|mps> and <mps|c+...c...|mps_aux> for symmetrised derivatives.
                        // The symmetrisation will be taken care for later
                        RDMEvaluator(mps_aux, mps);

                        // reset the current element back to 0
                        tst.data()[i](j,k) = 0.0;
                    }
                }
                else // one-site
                {
                    // do everything as above but operate on the MPSTensor directly rather than forming a two-site tensor
                    MPSTensor<Matrix, SymmGroup> & mpst_mod = mps_aux[site];
                    mpst_mod.make_left_paired();

                    mpst_mod.multiply_by_scalar(0.0);

                    for (int i = 0; i < mpst_mod.data().n_blocks(); i++)
                    for (int j = 0; j < mpst_mod.data()[i].num_rows(); j++)
                    for (int k = 0; k < mpst_mod.data()[i].num_cols(); k++)
                    {
                        this->ext_labels = { i, j, k };

                        mpst_mod.data()[i](j,k) = std::copysign(1.0, maquis::real(mps_aux[site].data()[i](j,k)));
                        mpst_mod.make_left_paired();
                        RDMEvaluator(mps_aux, mps);
                        mpst_mod.data()[i](j,k) = 0.0;
                    }
                }
            }
        private:
            int lr_site;

        };
    } // measurements

#endif