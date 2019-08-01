/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2019         Leon Freitag <lefreita@ethz.ch>
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
#ifndef NRDM_XVECTOR_H
#define NRDM_XVECTOR_H

// Measurement of transition density matrices for the variation of the MPS tensor
// expressed in nonredundant MPS parameters (XVector)
#include "dmrg/models/measurements/tagged_nrankrdm.h"
#include "dmrg/mp_tensors/xvector.h"
#include "dmrg/models/measurements/measurement_map.h"

    namespace measurements
    {
        template <class Matrix, class SymmGroup>
        class NRDMXVector : public TaggedNRankRDM<Matrix, SymmGroup>
        {
        typedef TaggedNRankRDM<Matrix, SymmGroup> base;
        typedef typename base::tag_vec tag_vec;
        typedef typename base::scaled_bond_term scaled_bond_term;
        typedef typename base::positions_type positions_type;
        typedef typename base::meas_type meas_type;

        using base::ext_labels;
        using base::do_measure_1rdm;
        using base::do_measure_2rdm;
        using base::lattice;

        public:
            // Specialization to call the correct constructor of the base class
            // TODO: Check if this compiles with disabled SU2U1/SU2U1PG!
            // symm_traits::HasSU2<SU2U1> and symm_traits::HasSU2<SU2U1PG> yields boost::true_type
            NRDMXVector(boost::true_type, std::string const& chkp, std::string const& aux_filename, std::string const& name_, const Lattice & lat,
                       boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_,
                       typename TermMakerSU2<Matrix, SymmGroup>::OperatorCollection const & op_collection_,
                       positions_type const& positions_ = positions_type(), bool twosite_ = true)
                       : base(name_, lat, tag_handler_, op_collection_, positions_), chkp_(chkp), aux_filename_(aux_filename) {};

            // 2U1 and other symmetry groups
            NRDMXVector(int lr_site_, boost::false_type, std::string const& chkp, std::string const& aux_filename, std::string const& name_, const Lattice & lat,
                       boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_,
                       tag_vec const & identities_, tag_vec const & fillings_, std::vector<scaled_bond_term> const& ops_,
                       bool half_only_, positions_type const& positions_ = positions_type(), bool twosite_ = true)
                       : base(name_, lat, tag_handler_, identities_, fillings_, ops_, half_only_, positions_), chkp_(chkp), aux_filename_(aux_filename) {};

            virtual void evaluate(MPS<Matrix, SymmGroup> const& ket_mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
            {
                lr::XVector<Matrix, SymmGroup> xvec;
                xvec.load(this->chkp_);  // Load the block matrix structure
                xvec.assign_mps(ket_mps);  // Specify the underlying MPS for the xvector

                if (!this->aux_filename_.empty()) // If we have an aux filename specified, load the data from the auxiliary text file (used in the MOLCAS interface)
                    xvec.update_from_textfile(this->aux_filename_);

                if (this->name() == "onerdmxvector")
                    measure_xvector_1rdm(ket_mps, xvec);

                if (this->name() == "twordmxvector")
                    measure_xvector_2rdm(ket_mps, xvec);

            }
        protected:
            measurement<Matrix, SymmGroup>* do_clone() const
            {
                return new NRDMXVector(*this);
            }

            // Measure transition densities wrt variation of the MPSTensors
            // The transition density wrt the variation of the MPSTensor is the sum of all transition densities
            // for each site, where the corresponding MPSTensor has been replaced by its variation
            void measure_xvector_1rdm(MPS<Matrix, SymmGroup> const& mps, lr::XVector<Matrix, SymmGroup> const & xvec)
            {
                MPS<Matrix, SymmGroup> mps_aux = mps;
                // Treat the first site separately to fill the measurement map
                assert(mps.length() > 0);
                mps_aux[0] = xvec.getB(0);
                onerdm_map<meas_type> rdm_full = do_measure_1rdm(mps_aux, mps);
                mps_aux[0] = mps[0];

                // Add the transition densities from other sites successively
                for (int i = 1; i < mps.length(); i++)
                {
                    // replace the MPSTensor for site i with variation
                    mps_aux[i] = xvec.getB(i);
                    // Measure the transition density for this site
                    onerdm_map<meas_type> rdm_partial = do_measure_1rdm(mps_aux, mps);
                    // Add the site transition density to the full density, element per element
                    for (typename onerdm_map<meas_type>::iterator it = rdm_partial.begin(); it != rdm_partial.end(); it++)
                        rdm_full.at(it->first) += it->second;
                    // Get back the original MPSTensor to prepare for the substitution of the next site
                    mps_aux[i] = mps[i];
                }

                // Convert the total map to the old-style results with labels and results
                std::tie(this->labels, this->vector_results) = measurements_details::convert_meas_map_to_old_results(rdm_full, lattice, ext_labels);
            }

            // and the same for 2-rdm (unfortunately copy-pasted from 1-rdm, maybe there's an easy way to make it more generic)
            void measure_xvector_2rdm(MPS<Matrix, SymmGroup> const& mps, lr::XVector<Matrix, SymmGroup> const & xvec)
            {
                MPS<Matrix, SymmGroup> mps_aux = mps;
                // Treat the first site separately to fill the measurement map
                assert(mps.length() > 0);
                mps_aux[0] = xvec.getB(0);
                twordm_map<meas_type> rdm_full = do_measure_2rdm(mps_aux, mps);
                mps_aux[0] = mps[0];

                // Add the transition densities from other sites successively
                for (int i = 1; i < mps.length(); i++)
                {
                    // replace the MPSTensor for site i with variation
                    mps_aux[i] = xvec.getB(i);
                    // Measure the transition density for this site
                    twordm_map<meas_type> rdm_partial = do_measure_2rdm(mps_aux, mps);
                    // Add the site transition density to the full density, element per element
                    for (typename twordm_map<meas_type>::iterator it = rdm_partial.begin(); it != rdm_partial.end(); it++)
                        rdm_full.at(it->first) += it->second;
                    // Get back the original MPSTensor to prepare for the substitution of the next site
                    mps_aux[i] = mps[i];
                }

                // Convert the total map to the old-style results with labels and results
                std::tie(this->labels, this->vector_results) = measurements_details::convert_meas_map_to_old_results(rdm_full, lattice, ext_labels);
            }

            std::string chkp_;
            std::string aux_filename_;
        };
    }
#endif