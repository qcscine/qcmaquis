/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MEASUREMENTS_LOCAL_H
#define MEASUREMENTS_LOCAL_H

#include <utility>

#include "dmrg/models/measurement.h"
#include "dmrg/models/meas_prepare.hpp"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/super_mpo.h"

namespace measurements {

    template <class Matrix, class SymmGroup>
    class local : public measurement<Matrix, SymmGroup> {
        using base = measurement<Matrix, SymmGroup>;
        using op_t = typename base::op_t;
        using generator = generate_mpo::MPOMaker<Matrix, SymmGroup>;
        using op_vec = std::vector<op_t>;
        using bond_element = std::vector<std::pair<op_vec, bool>>;
    public:

        local(std::string const& name_, Lattice  lat,
              op_vec const & identities_, op_vec const & fillings_,
              std::vector<bond_element> const& terms)
        : base(name_)
        , lattice(std::move(lat))
        , identities(identities_)
        , fillings(fillings_)
        , is_bond(true)
        , mpo_terms(terms)
        {
            this->cast_to_real = all_true(mpo_terms.begin(), mpo_terms.end(), static_cast<bool (*)(bond_element const&)>(&is_hermitian_meas));
        }

        local(std::string const& name_, Lattice  lat,
              op_vec const & identities_, op_vec const & fillings_,
              op_vec const& op)
        : base(name_)
        , lattice(std::move(lat))
        , identities(identities_)
        , fillings(fillings_)
        , is_bond(false)
        , site_term(op)
        , mpo_terms(std::vector<bond_element>(1, bond_element(1, std::make_pair(op, false)) ))
        {
            this->cast_to_real = is_hermitian_meas(site_term);
        }

        void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
        {
            this->vector_results.clear();
            this->labels.clear();
            this->labels_num.clear();

            using subcharge = typename SymmGroup::subcharge;
            if (!rmps || this->is_super_meas || is_bond) {
                evaluate_with_mpo(mps);
            } else {

                /// compute local reduced density matrices
                rmps.get().init();

                std::size_t L = mps.size();
                this->vector_results.reserve(this->vector_results.size() + L);
                this->labels.reserve(this->labels.size() + L);

                parallel::scheduler_balanced scheduler(L);

                for (typename Lattice::pos_t p = 0; p < L; ++p) {
                    parallel::guard proc(scheduler(p)); /// scheduling kernels

                    subcharge type = lattice.get_prop<subcharge>("type", p);
                    if (site_term[type].n_blocks() > 0) {

                        MPOTensor<Matrix, SymmGroup> temp;
                        temp.set(0, 0, site_term[type]);

                        MPSTensor<Matrix, SymmGroup> vec2
                        = contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(mps[p], rmps.get().left(p), rmps.get().right(p), temp);

                        typename MPS<Matrix, SymmGroup>::scalar_type res = mps[p].scalar_overlap(vec2);

                        this->vector_results.push_back(res);
                        this->labels.push_back( lattice.get_prop<std::string>("label", p) );
                        this->labels_num.push_back({p});
                    }
                }
            }
        }

    protected:
        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new local(*this);
        }

        void evaluate_with_mpo(MPS<Matrix, SymmGroup> const& mps)
        {
            using subcharge = typename SymmGroup::subcharge;
            using result_type = std::map<std::string, typename Matrix::value_type>;
            result_type res;

            typename MPS<Matrix, SymmGroup>::scalar_type nn;
            if (this->is_super_meas)
                nn = dm_trace(mps, this->phys_psi);

            /// collect results from all mpo terms, i.e. all requested combinations of operators.
            for (const auto& it : mpo_terms){
              using mpo_map = std::map<std::string, MPO<Matrix, SymmGroup>>;
              mpo_map mpos = meas_prepare::local<Matrix, SymmGroup>(lattice, identities, fillings, it);

              /// measure the value at each site / bond
              for (const auto& mit : mpos) {
                typename result_type::iterator match = res.find(mit.first);
                if (match == res.end())
                  boost::tie(match, boost::tuples::ignore) = res.insert( std::make_pair(mit.first, 0.) );

                if (!this->is_super_meas) {
                  match->second += (this->cast_to_real) ? maquis::real(expval(mps, mit.second)) : expval(mps, mit.second);
                } else {
                  MPS<Matrix, SymmGroup> super_mpo = mpo_to_smps(mit.second, this->phys_psi);
                  // static_cast needed for icpc 12.x
                  using overlap_func = typename MPS<Matrix, SymmGroup>::scalar_type (*)(const MPS<Matrix, SymmGroup> &, const MPS<Matrix, SymmGroup> &);
                  typename MPS<Matrix, SymmGroup>::scalar_type val = ::overlap(super_mpo, mps);
                  match->second += val/nn;
                }
              }

            }
            /// copy results to base
            this->vector_results.reserve(this->vector_results.size() + res.size());
            this->labels.reserve(this->labels.size() + res.size());
            for (const auto& it : res) {
                this->labels.push_back(it.first);
                this->vector_results.push_back(it.second);
            }
        }

    private:
        Lattice lattice;
        op_vec identities, fillings;
        bool is_bond;
        op_vec site_term;
        std::vector<bond_element> mpo_terms;
    };

}

#endif
