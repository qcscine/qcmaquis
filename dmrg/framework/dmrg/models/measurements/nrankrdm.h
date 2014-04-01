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
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/super_mpo.h"

#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>

namespace measurements {
    
    template <class Matrix, class SymmGroup>
    class NRankRDM : public measurement<Matrix, SymmGroup> {
        typedef measurement<Matrix, SymmGroup> base;
        typedef std::size_t size_type;
        typedef std::vector<block_matrix<Matrix, SymmGroup> > op_vec;
        typedef std::vector<size_type> positions_type;
        typedef Lattice::pos_t pos_t;
    
    public:
        NRankRDM(std::string const& name_, const Lattice & lat,
                 op_vec const & identities_, op_vec const & fillings_,
                 std::vector<std::pair<op_vec, bool> > const& ops_,
                 bool half_only_, bool nearest_neighbors_only,
                 positions_type const& positions_ = positions_type())
        : base(name_)
        , lattice(lat)
        , positions_first(positions_)
        , identities(identities_)
        , fillings(fillings_)
        , ops(ops_)
        , half_only(half_only_)
        , is_nn(nearest_neighbors_only)
        {
            pos_t extent = ops.size() > 2 ? lattice.size() : lattice.size()-1;
            if (positions_first.size() == 0)
                std::copy(boost::counting_iterator<int>(0), boost::counting_iterator<int>(extent),
                          back_inserter(positions_first));
            
            this->cast_to_real = is_hermitian_meas(ops);
        }
        
        void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
        {
            this->vector_results.clear();
            //this->labels.clear();

            if (ops.size() == 2)
                measure_correlation(mps, ops);
            else if (ops.size() == 4)
                measure_2rdm(mps, ops);
            else
                throw std::runtime_error("correlation measurements at the moment supported with 2 and 4 operators");
        }
        
    protected:
        typedef boost::shared_ptr<generate_mpo::CorrMakerBase<Matrix, SymmGroup> > maker_ptr;

        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new NRankRDM(*this);
        }
        
        void measure_correlation(MPS<Matrix, SymmGroup> const & mps,
                                 std::vector<std::pair<op_vec, bool> > const & ops,
                                 std::vector<size_type> const & order = std::vector<size_type>())
        {
            for (std::vector<std::size_t>::const_iterator it = positions_first.begin(); it != positions_first.end(); ++it) {
                #ifndef NDEBUG
                maquis::cout << "  site " << *it << std::endl;
                #endif
                
                maker_ptr dcorr(new generate_mpo::BgCorrMaker<Matrix, SymmGroup>(lattice, identities, fillings,
                                                                                 ops, std::vector<pos_t>(1, *it)));

                /// measure
                MPO<Matrix, SymmGroup> mpo = dcorr->create_mpo();
                std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct;
                dct = multi_expval(mps, mpo);
                
                /// save results and labels
                this->vector_results.reserve(this->vector_results.size() + dct.size());
                this->labels.reserve(this->labels.size() + dct.size());
                
                std::copy(dct.begin(), dct.end(), std::back_inserter(this->vector_results));
                
                //std::vector<std::vector<std::size_t> > num_labels = dcorr->numeric_labels();
                //std::vector<std::string> lbt = label_strings(lattice,  (order.size() > 0)
                //        ? detail::resort_labels(num_labels, order, is_nn) : num_labels );
                //std::copy(lbt.begin(), lbt.end(), std::back_inserter(this->labels));
            }
        }

        void measure_2rdm(MPS<Matrix, SymmGroup> const & mps,
                                 std::vector<std::pair<op_vec, bool> > const & ops,
                                 std::vector<size_type> const & order = std::vector<size_type>())
        {
            // TODO: test with ambient in due time
            #ifdef MAQUIS_OPENMP
            #pragma omp parallel for collapse(2)
            #endif
            for (size_t p1 = 0; p1 < lattice.size(); ++p1)
                for (size_t p2 = 0; p2 < lattice.size(); ++p2)
                {
                    size_t subref = std::min(p1, p2);
                    for (size_t p3 = subref; p3 < lattice.size(); ++p3)
                    { 
                        std::vector<pos_t> ref;
                        ref.push_back(p1); ref.push_back(p2); ref.push_back(p3);
                        maker_ptr dcorr(new generate_mpo::BgCorrMaker<Matrix, SymmGroup>(lattice, identities, fillings, ops, ref, true));

                        /// measure
                        MPO<Matrix, SymmGroup> mpo = dcorr->create_mpo();
                        std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct = multi_expval(mps, mpo);
                        
                        /// save results and labels
                        #ifdef MAQUIS_OPENMP
                        #pragma omp critical
                        #endif
                        {
                        this->vector_results.reserve(this->vector_results.size() + dct.size());
                        std::copy(dct.rbegin(), dct.rend(), std::back_inserter(this->vector_results));
                        }
                    }
                }
        }
        
    private:
        Lattice lattice;
        positions_type positions_first;
        op_vec identities, fillings;
        std::vector<std::pair<op_vec, bool> > ops;
        bool half_only, is_nn;
    };
}

#endif
