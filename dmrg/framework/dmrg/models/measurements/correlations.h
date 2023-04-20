/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MEASUREMENTS_CORRELATIONS_H
#define MEASUREMENTS_CORRELATIONS_H

#include "dmrg/models/measurement.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/super_mpo.h"

#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>

namespace measurements {

    namespace detail {
        typedef Lattice::pos_t pos_t;
        inline std::vector<std::vector<pos_t> >
        resort_labels (const std::vector<std::vector<pos_t> >& labels,
                       std::vector<pos_t> const & order,
                       bool is_nn=false)
        {
            std::vector<std::vector<pos_t> > ret(labels.size());
            for (int i=0; i<labels.size(); ++i) {
                #ifndef NDEBUG
                if (is_nn) assert(2*order.size() == labels[i].size());
                else       assert(order.size() == labels[i].size());
                #endif
                ret[i].resize(labels[i].size());

                for (int j=0; j<order.size(); ++j) {
                    if (is_nn) {
                        ret[i][order[j]] = labels[i][2*j];
                        ret[i][order[j]+1] = labels[i][2*j+1];
                    } else {
                        ret[i][order[j]] = labels[i][j];
                    }
                }
            }
            return ret;
        }


        template <class Matrix, class SymmGroup>
        class CorrPermutator
        {
        public:
            typedef Lattice::pos_t pos_t;
            typedef std::size_t size_t;
            typedef typename measurement<Matrix, SymmGroup>::op_t op_t;
            typedef std::pair<std::vector<op_t>, bool> inner_t;
            typedef std::vector<inner_t> value_t;
            CorrPermutator (value_t const & ops, bool is_nn)
            {
                std::vector<pos_t> ord;
                if (is_nn)
                    for (pos_t i=0; i<ops.size(); i+=2) ord.push_back(i);
                else
                    for (pos_t i=0; i<ops.size(); ++i) ord.push_back(i);

                do {
                    std::vector<pos_t> check_pre;
                    for (pos_t i=0; i<ops.size(); ++i) {
                        if (ops[i].second) {
                            if (!is_nn)
                                check_pre.push_back(ord[i]);
                            else
                                if (i % 2 == 0)
                                    check_pre.push_back(ord[i]);
                                else
                                    check_pre.push_back(ord[i]+1);
                        }
                    }

                    cmp_with_prefactor cmp;
                    cmp.prefactor = 1.;
                    std::sort(check_pre.begin(), check_pre.end(), cmp);

                    value_t tmp;
                    for (pos_t i=0; i<ord.size(); ++i) {
                        tmp.push_back( ops[ord[i]] );
                        if (is_nn) {
                            tmp.push_back( ops[ord[i]+1] );
                        }
                    }

                    for (int type=0; type<tmp[0].first.size(); ++type)
                        tmp[0].first[type] *= cmp.prefactor;

                    perm.push_back(tmp);
                    orders.push_back(ord);
                } while (std::next_permutation(ord.begin(), ord.end()));
            }

            size_t size() const {return perm.size();}
            value_t operator[](size_t i) const {return perm[i];}
            std::vector<pos_t> order(size_t i) const {return orders[i];}

        private:
            std::vector<value_t> perm;
            std::vector<std::vector<pos_t> > orders;
        };
    } /// namespace detail



    template <class Matrix, class SymmGroup>
    class correlations : public measurement<Matrix, SymmGroup> {
        typedef measurement<Matrix, SymmGroup> base;
        typedef typename base::op_t op_t;
        typedef Lattice::pos_t pos_t;
        typedef std::vector<op_t> op_vec;
        typedef std::vector<pos_t> positions_type;

    public:
        correlations(std::string const& name_, const Lattice & lat,
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
            if (positions_first.size() == 0)
                std::copy(boost::counting_iterator<int>(0), boost::counting_iterator<int>(lattice.size()-(ops.size()-1)),
                          back_inserter(positions_first));

            this->cast_to_real = is_hermitian_meas(ops);
        }

        void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
        {
            this->vector_results.clear();
            this->labels.clear();
            this->labels_num.clear();

            if (half_only) {
                measure_correlation(mps, ops);
            } else {
                detail::CorrPermutator<Matrix, SymmGroup> perm(ops, is_nn);
                for (int i=0; i<perm.size(); ++i) {
                    measure_correlation(mps, perm[i], perm.order(i));
                }
            }
        }

    protected:
        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new correlations(*this);
        }

        void measure_correlation(MPS<Matrix, SymmGroup> const & mps,
                                 std::vector<std::pair<op_vec, bool> > const & ops,
                                 std::vector<pos_t> const & order = std::vector<pos_t>())
        {
            typedef std::shared_ptr<generate_mpo::CorrMakerBase<Matrix, SymmGroup> > maker_ptr;

            for (std::vector<pos_t>::const_iterator it = positions_first.begin(); it != positions_first.end(); ++it) {
                if (*it >= lattice.size()-(ops.size()-1))
                    throw std::runtime_error("cannot measure correlation with first operator at p="+boost::lexical_cast<std::string>(*it)+".");
                #ifndef NDEBUG
                maquis::cout << "  site " << *it << std::endl;
                #endif

                /// initialize correct maker
                maker_ptr dcorr;
                if (is_nn) dcorr.reset(new generate_mpo::CorrMakerNN<Matrix, SymmGroup>(lattice, identities, fillings, ops, *it) );
                else       dcorr.reset(new generate_mpo::CorrMaker<Matrix, SymmGroup>(lattice, identities, fillings, ops, *it) );

                /// measure
                MPO<Matrix, SymmGroup> mpo = dcorr->create_mpo();
                std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct;
                if (!this->is_super_meas)
                    dct = multi_expval(mps, mpo);
                else {
                    typename MPS<Matrix, SymmGroup>::scalar_type nn = dm_trace(mps, this->phys_psi);
                    MPS<Matrix, SymmGroup> super_mpo = mpo_to_smps(mpo, this->phys_psi);
                    dct = multi_overlap(super_mpo, mps);
                    for (int i=0; i<dct.size(); ++i)
                        dct[i] /= nn;
                }

                /// save results and labels
                this->vector_results.reserve(this->vector_results.size() + dct.size());
                this->labels.reserve(this->labels.size() + dct.size());
                this->labels_num.reserve(this->labels_num.size() + dct.size());

                std::copy(dct.begin(), dct.end(), std::back_inserter(this->vector_results));

                std::vector<std::vector<pos_t> > num_labels = (order.size() > 0) ? detail::resort_labels(dcorr->numeric_labels(), order, is_nn) : dcorr->numeric_labels();

                std::copy(num_labels.begin(), num_labels.end(), std::back_inserter(this->labels_num));

                std::vector<std::string> lbt = label_strings(lattice, num_labels );
                std::copy(lbt.begin(), lbt.end(), std::back_inserter(this->labels));
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
