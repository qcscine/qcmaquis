/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
 *               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef MEASUREMENTS_NRANKRDM_H
#define MEASUREMENTS_NRANKRDM_H

#include "dmrg/models/measurement.h"
#include "dmrg/models/measurements/correlations.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/super_mpo.h"

#include <algorithm>
#include <functional>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/filesystem.hpp>

namespace measurements {
    
    template <class Matrix, class SymmGroup>
    class NRankRDM : public measurement<Matrix, SymmGroup> {
        typedef measurement<Matrix, SymmGroup> base;
        typedef std::vector<block_matrix<Matrix, SymmGroup> > op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
        typedef Lattice::pos_t pos_t;
        typedef std::vector<pos_t> positions_type;
    
    public:
        NRankRDM(std::string const& name_, const Lattice & lat,
                 op_vec const & identities_, op_vec const & fillings_,
                 std::vector<bond_element> const& ops_,
                 bool half_only_, bool nearest_neighbors_only,
                 positions_type const& positions_ = positions_type(),
                 std::string const& ckp_ = std::string(""))
        : base(name_)
        , lattice(lat)
        , positions_first(positions_)
        , identities(identities_)
        , fillings(fillings_)
        , ops(ops_)
        , half_only(half_only_)
        , is_nn(nearest_neighbors_only)
        , bra_ckp(ckp_)
        {
            pos_t extent = ops.size() > 2 ? lattice.size() : lattice.size()-1;
            if (positions_first.size() == 0)
                std::copy(boost::counting_iterator<pos_t>(0), boost::counting_iterator<pos_t>(extent),
                          back_inserter(positions_first));
            
            this->cast_to_real = is_hermitian_meas(ops[0]);
        }
        
        void evaluate(MPS<Matrix, SymmGroup> const& ket_mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
        {
            this->vector_results.clear();
            this->labels.clear();

            MPS<Matrix, SymmGroup> bra_mps;
            if (bra_ckp != "") {
                if(boost::filesystem::exists(bra_ckp))
                    load(bra_ckp, bra_mps);
                else
                    throw std::runtime_error("The bra checkpoint file " + bra_ckp + " was not found\n");
            }

            if (ops[0].size() == 2)
                measure_correlation(bra_mps, ket_mps, ops);
            else if (ops[0].size() == 4)
                measure_2rdm(bra_mps, ket_mps, ops);
            else if (ops[0].size() == 6)
                measure_3rdm(bra_mps, ket_mps, ops);
            else
                throw std::runtime_error("correlation measurements at the moment supported with up to 6 operators, you have "
                                          + boost::lexical_cast<std::string>(ops[0].size()) + "\n");
        }
        
    protected:
        typedef boost::shared_ptr<generate_mpo::CorrMakerBase<Matrix, SymmGroup> > maker_ptr;

        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new NRankRDM(*this);
        }
        
        void measure_correlation(MPS<Matrix, SymmGroup> const & dummy_bra_mps,
                                 MPS<Matrix, SymmGroup> const & ket_mps,
                                 std::vector<bond_element> const & ops,
                                 std::vector<pos_t> const & order = std::vector<pos_t>())
        {
            // Test if a separate bra state has been specified
            bool bra_neq_ket = (dummy_bra_mps.length() > 0);
            MPS<Matrix, SymmGroup> const & bra_mps = (bra_neq_ket) ? dummy_bra_mps : ket_mps;

            // TODO: test with ambient in due time
            #ifdef MAQUIS_OPENMP
            #pragma omp parallel for
            #endif
            for (std::size_t i = 0; i < positions_first.size(); ++i) {
                pos_t p = positions_first[i];
                #ifndef NDEBUG
                maquis::cout << "  site " << p << std::endl;
                #endif
                
                maker_ptr dcorr(new generate_mpo::BgCorrMaker<Matrix, SymmGroup>(lattice, identities, fillings,
                                                                                 ops[0], std::vector<pos_t>(1, p)));
                // measure
                MPO<Matrix, SymmGroup> mpo = dcorr->create_mpo();
                std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct = multi_expval(bra_mps, ket_mps, mpo);
                
                std::vector<std::vector<pos_t> > num_labels = dcorr->numeric_labels();
                std::vector<std::string> lbt = label_strings(lattice,  (order.size() > 0)
                                            ? detail::resort_labels(num_labels, order, is_nn) : num_labels );
                // save results and labels
                #ifdef MAQUIS_OPENMP
                #pragma omp critical
                #endif
                {
                this->vector_results.reserve(this->vector_results.size() + dct.size());
                std::copy(dct.begin(), dct.end(), std::back_inserter(this->vector_results));

                this->labels.reserve(this->labels.size() + dct.size());
                std::copy(lbt.begin(), lbt.end(), std::back_inserter(this->labels));
                }
            }
        }

        void measure_2rdm(MPS<Matrix, SymmGroup> const & dummy_bra_mps,
                          MPS<Matrix, SymmGroup> const & ket_mps,
                          std::vector<bond_element> const & ops,
                          std::vector<pos_t> const & order = std::vector<pos_t>())
        {
            // Test if a separate bra state has been specified
            bool bra_neq_ket = (dummy_bra_mps.length() > 0);
            MPS<Matrix, SymmGroup> const & bra_mps = (bra_neq_ket) ? dummy_bra_mps : ket_mps;

            // TODO: test with ambient in due time
            #ifdef MAQUIS_OPENMP
            #pragma omp parallel for collapse(2)
            #endif
            for (pos_t p1 = 0; p1 < lattice.size(); ++p1)
            for (pos_t p2 = 0; p2 < lattice.size(); ++p2)
            {
                // Permutation symmetry for bra == ket: ijkl == jilk == klji == lkji
                pos_t subref = std::min(p1, p2);

                // if bra != ket, pertmutation symmetry is only ijkl == jilk
                if (bra_neq_ket)
                    pos_t subref = 0;

                for (pos_t p3 = subref; p3 < lattice.size(); ++p3)
                { 
                    // Measurement positions p1,p2,p3 are fixed, p4 is handled by the MPO (synmpo)
                    std::vector<pos_t> ref;
                    ref.push_back(p1); ref.push_back(p2); ref.push_back(p3);

                    maker_ptr dcorr(new generate_mpo::BgCorrMaker<Matrix, SymmGroup>(lattice, identities, fillings, ops[0], ref, true));
                    MPO<Matrix, SymmGroup> mpo = dcorr->create_mpo();
                    std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct = multi_expval(bra_mps, ket_mps, mpo);

                    // Loop over operator terms that are measured synchronously and added together
                    // Used e.g. for the four spin combos of the 2-RDM
                    for (std::size_t synop = 1; synop < ops.size(); ++synop) {
                        maker_ptr syndcorr(new generate_mpo::BgCorrMaker<Matrix, SymmGroup>(lattice, identities, fillings, ops[synop], ref, true));

                        // measure
                        MPO<Matrix, SymmGroup> synmpo = syndcorr->create_mpo();
                        std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> syndct = multi_expval(bra_mps, ket_mps, synmpo);

                        // add synchronous terms
                        std::transform(syndct.begin(), syndct.end(), dct.begin(), dct.begin(),
                                       std::plus<typename MPS<Matrix, SymmGroup>::scalar_type>());
                    }
                    
                    std::vector<std::vector<pos_t> > num_labels = dcorr->numeric_labels();
                    std::vector<std::string> lbt = label_strings(lattice,  (order.size() > 0)
                                                ? detail::resort_labels(num_labels, order, is_nn) : num_labels );
                    // save results and labels
                    #ifdef MAQUIS_OPENMP
                    #pragma omp critical
                    #endif
                    {
                    this->vector_results.reserve(this->vector_results.size() + dct.size());
                    std::copy(dct.rbegin(), dct.rend(), std::back_inserter(this->vector_results));

                    this->labels.reserve(this->labels.size() + dct.size());
                    std::copy(lbt.rbegin(), lbt.rend(), std::back_inserter(this->labels));
                    }
                }
            }
        }

        void measure_3rdm(MPS<Matrix, SymmGroup> const & dummy_bra_mps,
                          MPS<Matrix, SymmGroup> const & ket_mps,
                          std::vector<bond_element> const & ops,
                          std::vector<pos_t> const & order = std::vector<pos_t>())
        {
            // Test if a separate bra state has been specified
            bool bra_neq_ket = (dummy_bra_mps.length() > 0);
            MPS<Matrix, SymmGroup> const & bra_mps = (bra_neq_ket) ? dummy_bra_mps : ket_mps;

            // TODO: test with ambient in due time
            #ifdef MAQUIS_OPENMP
            #pragma omp parallel for collapse(2)
            #endif
            //throw std::runtime_error("stefan: 3RDM measurement not yet supported");
            for (pos_t p1 = 0; p1 < lattice.size(); ++p1)
            for (pos_t p2 = 0; p2 < lattice.size(); ++p2)
            for (pos_t p3 = 0; p3 < lattice.size(); ++p3)
            for (pos_t p4 = 0; p4 < lattice.size(); ++p4)
            for (pos_t p5 = 0; p5 < lattice.size(); ++p5)
            {
                // Measurement positions p1,p2,p3,p4,p5 are fixed
                std::vector<pos_t> ref;
                ref.push_back(p1); ref.push_back(p2); ref.push_back(p3); ref.push_back(p4); ref.push_back(p5);

                // ops[0] is the first set of 6 operators that will be measured 
                maker_ptr dcorr(new generate_mpo::BgCorrMaker<Matrix, SymmGroup>(lattice, identities, fillings, ops[0], ref, true));
                MPO<Matrix, SymmGroup> mpo = dcorr->create_mpo();
                std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct = multi_expval(bra_mps, ket_mps, mpo);

                // Loop over operator terms that are measured synchronously and added together
                // Used e.g. for the four spin combos of the 2-RDM
                for (std::size_t synop = 1; synop < ops.size(); ++synop) {
                    maker_ptr syndcorr(new generate_mpo::BgCorrMaker<Matrix, SymmGroup>(lattice, identities, fillings, ops[synop], ref, true));

                    // measure
                    MPO<Matrix, SymmGroup> synmpo = syndcorr->create_mpo();
                    std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> syndct = multi_expval(bra_mps, ket_mps, synmpo);

                    // add synchronous terms
                    std::transform(syndct.begin(), syndct.end(), dct.begin(), dct.begin(),
                                   std::plus<typename MPS<Matrix, SymmGroup>::scalar_type>());
                }
                
                // the label consists of p1,...,p5; the 6th label is implicit, because terms are in order
                std::vector<std::vector<pos_t> > num_labels = dcorr->numeric_labels();
                std::vector<std::string> lbt = label_strings(lattice,  (order.size() > 0)
                                            ? detail::resort_labels(num_labels, order, is_nn) : num_labels );
                // save results and labels
                #ifdef MAQUIS_OPENMP
                #pragma omp critical
                #endif
                {
                this->vector_results.reserve(this->vector_results.size() + dct.size());
                std::copy(dct.rbegin(), dct.rend(), std::back_inserter(this->vector_results));

                this->labels.reserve(this->labels.size() + dct.size());
                std::copy(lbt.rbegin(), lbt.rend(), std::back_inserter(this->labels));
                }
            }
        }
        
    private:
        Lattice lattice;
        positions_type positions_first;
        op_vec identities, fillings;
        std::vector<bond_element> ops;
        bool half_only, is_nn;

        std::string bra_ckp;
    };
}

#endif
