/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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

#ifndef MEASUREMENTS_TAGGED_NRANKRDM_H
#define MEASUREMENTS_TAGGED_NRANKRDM_H

#include <algorithm>
#include <functional>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/filesystem.hpp>

#include "dmrg/models/term_descriptor.h"
#include "dmrg/models/measurement.h"

namespace measurements {
    
    template <class Matrix, class SymmGroup>
    class TaggedNRankRDM : public measurement<Matrix, SymmGroup> {

        typedef measurement<Matrix, SymmGroup> base;

        typedef typename Model<Matrix, SymmGroup>::term_descriptor term_descriptor;

        typedef Lattice::pos_t pos_t;
        typedef std::vector<pos_t> positions_type;

        typedef typename base::op_t op_t;
        typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
        typedef std::vector<tag_type> tag_vec;
        typedef std::vector<tag_vec> bond_term;
    
    public:
        TaggedNRankRDM(std::string const& name_, const Lattice & lat,
                       boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_,
                       tag_vec const & identities_, tag_vec const & fillings_, std::vector<bond_term> const& ops_,
                       bool half_only_, positions_type const& positions_ = positions_type(),
                       std::string const& ckp_ = std::string(""))
        : base(name_)
        , lattice(lat)
        , tag_handler(tag_handler_)
        , positions_first(positions_)
        , identities(identities_)
        , fillings(fillings_)
        , operator_terms(ops_)
        , bra_ckp(ckp_)
        {
            pos_t extent = operator_terms.size() > 2 ? lattice.size() : lattice.size()-1;
            if (positions_first.size() == 0)
                std::copy(boost::counting_iterator<pos_t>(0), boost::counting_iterator<pos_t>(extent),
                          back_inserter(positions_first));
            
            //this->cast_to_real = is_hermitian_meas(ops[0]);
            this->cast_to_real = false;
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

            if (operator_terms[0].size() == 2)
                measure_correlation(bra_mps, ket_mps);
            else if (operator_terms[0].size() == 4)
                measure_2rdm(bra_mps, ket_mps);
            else if (operator_terms[0].size() == 6)
                measure_3rdm(bra_mps, ket_mps);
            else
                throw std::runtime_error("correlation measurements at the moment supported with 2, 4 and 6 operators, size is "
                                          + boost::lexical_cast<std::string>(operator_terms[0].size()));
        }
        
    protected:

        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new TaggedNRankRDM(*this);
        }
        
        void measure_correlation(MPS<Matrix, SymmGroup> const & dummy_bra_mps,
                                 MPS<Matrix, SymmGroup> const & ket_mps,
                                 std::vector<pos_t> const & order = std::vector<pos_t>())
        {
            // Test if a separate bra state has been specified
            bool bra_neq_ket = (dummy_bra_mps.length() > 0);
            MPS<Matrix, SymmGroup> const & bra_mps = (bra_neq_ket) ? dummy_bra_mps : ket_mps;

            #ifdef MAQUIS_OPENMP
            #pragma omp parallel for
            #endif
            for (std::size_t i = 0; i < positions_first.size(); ++i) {
                pos_t p1 = positions_first[i];
                boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_local(new TagHandler<Matrix, SymmGroup>(*tag_handler));

                std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct;
                std::vector<std::vector<pos_t> > num_labels;
                for (pos_t p2 = p1+1; p2 < lattice.size(); ++p2)
                { 
                    pos_t pos_[2] = {p1, p2};
                    std::vector<pos_t> positions(pos_, pos_ + 2);

                    tag_vec operators(2);
                    operators[0] = operator_terms[0][0][lattice.get_prop<typename SymmGroup::subcharge>("type", p1)];
                    operators[1] = operator_terms[0][1][lattice.get_prop<typename SymmGroup::subcharge>("type", p2)];
                    
                    //term_descriptor term = generate_mpo::arrange_operators(tag_handler, positions, operators);
                    MPO<Matrix, SymmGroup> mpo = generate_mpo::make_1D_mpo(positions, operators, identities, fillings, tag_handler_local, lattice);
                    typename MPS<Matrix, SymmGroup>::scalar_type value = expval(bra_mps, ket_mps, mpo);

                    dct.push_back(value);
                    num_labels.push_back(positions);
                }

                std::vector<std::string> lbt = label_strings(lattice,  (order.size() > 0)
                                            ? detail::resort_labels(num_labels, order, false) : num_labels );
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

        void measure_2rdm(MPS<Matrix, SymmGroup> const & dummy_bra_mps,
                          MPS<Matrix, SymmGroup> const & ket_mps,
                          std::vector<pos_t> const & order = std::vector<pos_t>())
        {
            // Test if a separate bra state has been specified
            bool bra_neq_ket = (dummy_bra_mps.length() > 0);
            MPS<Matrix, SymmGroup> const & bra_mps = (bra_neq_ket) ? dummy_bra_mps : ket_mps;

            #ifdef MAQUIS_OPENMP
            #pragma omp parallel for collapse(1)
            #endif
            for (pos_t p1 = 0; p1 < lattice.size(); ++p1)
            for (pos_t p2 = 0; p2 < lattice.size(); ++p2)
            {
                boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_local(new TagHandler<Matrix, SymmGroup>(*tag_handler));

                // Permutation symmetry for bra == ket: pqrs == rspq == qpsr == srqp
                pos_t subref = std::min(p1, p2);

                // if bra != ket, pertmutation symmetry is only pqrs == qpsr
                if (bra_neq_ket)
                    pos_t subref = 0;

                for (pos_t p3 = subref; p3 < lattice.size(); ++p3)
                { 
                    std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct;
                    std::vector<std::vector<pos_t> > num_labels;

                    for (pos_t p4 = p3; p4 < lattice.size(); ++p4)
                    { 
                        pos_t pos_[4] = {p1, p2, p3, p4};
                        std::vector<pos_t> positions(pos_, pos_ + 4);

                        // Loop over operator terms that are measured synchronously and added together
                        // Used e.g. for the four spin combos of the 2-RDM
                        typename MPS<Matrix, SymmGroup>::scalar_type value = 0;
                        for (std::size_t synop = 0; synop < operator_terms.size(); ++synop) {

                            tag_vec operators(4);
                            operators[0] = operator_terms[synop][0][lattice.get_prop<typename SymmGroup::subcharge>("type", p1)];
                            operators[1] = operator_terms[synop][1][lattice.get_prop<typename SymmGroup::subcharge>("type", p2)];
                            operators[2] = operator_terms[synop][2][lattice.get_prop<typename SymmGroup::subcharge>("type", p3)];
                            operators[3] = operator_terms[synop][3][lattice.get_prop<typename SymmGroup::subcharge>("type", p4)];
                            
                            //term_descriptor term = generate_mpo::arrange_operators(tag_handler, positions, operators);
                            MPO<Matrix, SymmGroup> mpo = generate_mpo::make_1D_mpo(positions, operators, identities, fillings, tag_handler_local, lattice);
                            value += expval(bra_mps, ket_mps, mpo);
                        }

                        dct.push_back(value);
                        num_labels.push_back(positions);
                    }

                    std::vector<std::string> lbt = label_strings(lattice,  (order.size() > 0)
                                                ? detail::resort_labels(num_labels, order, false) : num_labels );

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
                          std::vector<pos_t> const & order = std::vector<pos_t>())
        {
            // Test if a separate bra state has been specified bool bra_neq_ket = (dummy_bra_mps.length() > 0);
            bool bra_neq_ket = (dummy_bra_mps.length() > 0);
            MPS<Matrix, SymmGroup> const & bra_mps = (bra_neq_ket) ? dummy_bra_mps : ket_mps;

            #ifdef MAQUIS_OPENMP
            #pragma omp parallel for collapse(3)
            #endif
            for (pos_t p1 = 0; p1 < lattice.size(); ++p1)
            for (pos_t p2 = 0; p2 < lattice.size(); ++p2)
            for (pos_t p3 = 0; p3 < lattice.size(); ++p3)
            for (pos_t p4 = 0; p4 < lattice.size(); ++p4)
            {
                boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_local(new TagHandler<Matrix, SymmGroup>(*tag_handler));

                for (pos_t p5 = 0; p5 < lattice.size(); ++p5)
                { 
                    std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct;
                    std::vector<std::vector<pos_t> > num_labels;

                    for (pos_t p6 = 0; p6 < lattice.size(); ++p6)
                    { 
                        pos_t pos_[6] = {p1, p2, p3, p4, p5, p6};
                        std::vector<pos_t> positions(pos_, pos_ + 6);

                        // Loop over operator terms that are measured synchronously and added together
                        // Used e.g. for the four spin combos of the 2-RDM
                        typename MPS<Matrix, SymmGroup>::scalar_type value = 0;
                        for (std::size_t synop = 0; synop < operator_terms.size(); ++synop) {

                            tag_vec operators(6);
                            operators[0] = operator_terms[synop][0][lattice.get_prop<typename SymmGroup::subcharge>("type", p1)];
                            operators[1] = operator_terms[synop][1][lattice.get_prop<typename SymmGroup::subcharge>("type", p2)];
                            operators[2] = operator_terms[synop][2][lattice.get_prop<typename SymmGroup::subcharge>("type", p3)];
                            operators[3] = operator_terms[synop][3][lattice.get_prop<typename SymmGroup::subcharge>("type", p4)];
                            operators[4] = operator_terms[synop][4][lattice.get_prop<typename SymmGroup::subcharge>("type", p5)];
                            operators[5] = operator_terms[synop][5][lattice.get_prop<typename SymmGroup::subcharge>("type", p6)];
                            
                            //term_descriptor term = generate_mpo::arrange_operators(tag_handler, positions, operators);
                            MPO<Matrix, SymmGroup> mpo = generate_mpo::make_1D_mpo(positions, operators, identities, fillings, tag_handler_local, lattice);
                            value += expval(bra_mps, ket_mps, mpo);
                        }

                        dct.push_back(value);
                        num_labels.push_back(positions);
                    }

                    std::vector<std::string> lbt = label_strings(lattice,  (order.size() > 0)
                                                ? detail::resort_labels(num_labels, order, false) : num_labels );

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
        
    private:
        Lattice lattice;
        boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
        positions_type positions_first;
        tag_vec identities, fillings;
        std::vector<bond_term> operator_terms;

        std::string bra_ckp;
    };
}

#endif
