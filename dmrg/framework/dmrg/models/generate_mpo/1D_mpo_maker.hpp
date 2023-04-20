/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef ONED_MPO_MAKER_H
#define ONED_MPO_MAKER_H

#include <utility>
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/models/model.h"
#include "dmrg/models/term_descriptor.h"
#include "dmrg/models/generate_mpo/utils.hpp"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/OperatorHandlers/OpTable.h"

namespace generate_mpo
{
    // currently not used
    template <class Matrix, class SymmGroup>
    std::vector<std::pair<Lattice::pos_t, std::pair<std::vector<typename OPTable<Matrix, SymmGroup>::op_t>, bool> > >
    arrange_operators(std::vector<Lattice::pos_t> const & positions,
                      std::vector<std::pair<std::vector<typename OPTable<Matrix, SymmGroup>::op_t>, bool> > const & operators)

    {
        // input: list of positions and operators
        // output: list of (position, operator)-pairs, sorted, unique positions with operators multiplied
        typedef Lattice::pos_t pos_t;
        typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;
        typedef std::pair<std::vector<op_t>, bool> site_ops_t;
        typedef std::pair<pos_t, site_ops_t> pos_op_t;

        std::vector<pos_op_t> pos_ops;
        // arrange position / operators in pairs
        std::transform(positions.begin(), positions.end()-1, operators.begin(), std::back_inserter(pos_ops),
        std::make_pair<pos_t const&, site_ops_t const&>);
                       // boost::bind(static_cast<pos_op_t(*)(pos_t const&, site_ops_t const&)>
                       // (std::make_pair<pos_t, site_ops_t>), boost::lambda::_1, boost::lambda::_2));

        std::stable_sort(pos_ops.begin(), pos_ops.end(), compare<pos_op_t>);

        for (size_t opnr = 0; opnr < pos_ops.size(); )
        {
            site_ops_t product = pos_ops[opnr].second;
            size_t range_end = opnr+1;

            // while the next operator is still on the same site
            while (range_end < pos_ops.size() && pos_ops[range_end].first == pos_ops[opnr].first) {
                // multiply operators for all irreps (types)
                for (size_t type=0; type < pos_ops[opnr].second.first.size(); ++type) {
                    op_t tmp;
                    gemm(pos_ops[range_end].second.first[type], product.first[type], tmp);
                    product.first[type] = tmp;
                }
                // set the fermion or boson for the product operator
                product.second = (pos_ops[range_end].second.second != product.second);
                range_end++;
            }
            for (size_t r=opnr; r < range_end; ++r)
                pos_ops[r].second = product;

            opnr = range_end;
        }
        return pos_ops;
    }


    template <class Matrix, class SymmGroup>
    term_descriptor<typename Matrix::value_type>
    arrange_operators(std::vector<Lattice::pos_t> const & positions,
                      std::vector<typename OPTable<Matrix, SymmGroup>::tag_type> const & operators,
                      std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler)
    {
        // input: list of positions and operators
        // output: list of (position, operator)-pairs, sorted, unique positions with operators multiplied

        assert(positions.size() == operators.size());

        typedef Lattice::pos_t pos_t;
        typedef typename Matrix::value_type value_type;
        typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;

        typedef term_descriptor<value_type> term_descriptor;

        typedef std::pair<pos_t, tag_type> pos_op_t;
        std::vector<pos_op_t> pos_ops;

        // determine the phase
        int phase = 1;
        {
            int inv_count=0, n=positions.size();
            for(pos_t c1 = 0; c1 < n - 1; c1++)
                for(pos_t c2 = c1+1; c2 < n; c2++)
                    if(positions[c1] > positions[c2]
                       && tag_handler->is_fermionic(operators[c1])
                       && tag_handler->is_fermionic(operators[c2])) inv_count++;

            if (inv_count % 2 == 1) phase = -1;
        }

        // arrange position / operators in pairs
        std::transform(positions.begin(), positions.end(), operators.begin(), std::back_inserter(pos_ops),
                       std::make_pair<pos_t const&, tag_type const&>);
        std::stable_sort(pos_ops.begin(), pos_ops.end(), generate_mpo::compare<pos_op_t>);

        term_descriptor term;
        term.coeff = phase;

        for (size_t opnr = 0; opnr < pos_ops.size(); )
        {
            tag_type product = pos_ops[opnr].second;
            size_t range_end = opnr + 1;

            // while the next operator is still on the same site
            while (range_end < pos_ops.size() && pos_ops[range_end].first == pos_ops[opnr].first) {
                value_type scale = 1.;
                boost::tie(product, scale) = tag_handler->get_product_tag(pos_ops[range_end].second, product);
                term.coeff *= scale;
                range_end++;
            }

            term.push_back( std::make_pair(pos_ops[opnr].first, product) );

            opnr = range_end;
        }
        return term;
    }

    namespace detail {
        template <int N, typename Tuple>
        inline typename boost::tuples::access_traits<typename boost::tuples::element<N, Tuple>::type>::const_type
        get(Tuple const& t)
        {
             return boost::tuples::get<N>(t) ;
        }
    }

    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup>
    sign_and_fill(term_descriptor<typename Matrix::value_type> term,
                  std::vector<typename OPTable<Matrix, SymmGroup>::tag_type> const & ident,
                  std::vector<typename OPTable<Matrix, SymmGroup>::tag_type> const & fill,
                  std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler,
                  Lattice const & lat)
    {
        // after arrange operators, expand term to the full site-list

        typedef typename SymmGroup::subcharge sc;
        typedef Lattice::pos_t pos_t;
        typedef term_descriptor<typename Matrix::value_type> term_descriptor;
        typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;

        MPO<Matrix, SymmGroup> ret(lat.size());

        bool carry_sign = false;
        for (pos_t p = 0; p < lat.size(); ++p)
        {
            typename MPOTensor<Matrix, SymmGroup>::prempo_t prempo;

            // check if there is a non-trivial operator on site p
            auto it = std::find_if(term.begin(), term.end(), [&p](const auto& iInput) {
                return iInput.first == p;
            });

            if (it != term.end()) // if yes
            {
                // multiply with fill if necessary
                tag_type product = get<1>(*it);
                if (tag_handler->is_fermionic(product) != carry_sign)
                {
                    typename Matrix::value_type scale;
                    boost::tie(product, scale) = tag_handler->get_product_tag(fill[lat.get_prop<sc>("type", p)], product);
                    term.coeff *= scale;
                }

                // update the phase
                if (tag_handler->is_fermionic(product)) carry_sign = !carry_sign;

                prempo.push_back(boost::make_tuple(0,0, product, 1.));
                ret[p] = MPOTensor<Matrix, SymmGroup>(1, 1, prempo, tag_handler->get_operator_table());
            }

            else if (carry_sign) // if no
            {
                prempo.push_back(boost::make_tuple(0,0, fill[lat.get_prop<sc>("type", p)], 1.));
                ret[p] = MPOTensor<Matrix, SymmGroup>(1, 1, prempo, tag_handler->get_operator_table());
            }

            else
            {
                prempo.push_back(boost::make_tuple(0,0, ident[lat.get_prop<sc>("type", p)], 1.));
                ret[p] = MPOTensor<Matrix, SymmGroup>(1, 1, prempo, tag_handler->get_operator_table());
            }

        }

        // put the term scale on the first non-trivial operator
        ret[term.position(0)].multiply_by_scalar(term.coeff);
        return ret;
    }

    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup>
    make_1D_mpo(std::vector<Lattice::pos_t> const & positions,
                std::vector<typename OPTable<Matrix, SymmGroup>::tag_type> const & operators,
                std::vector<typename OPTable<Matrix, SymmGroup>::tag_type> const & ident,
                std::vector<typename OPTable<Matrix, SymmGroup>::tag_type> const & fill,
                std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler,
                Lattice const & lat, typename Matrix::value_type scale = 1)
    {
        term_descriptor<typename Matrix::value_type> term = arrange_operators(positions, operators, tag_handler);
        term.coeff *= scale;
        return sign_and_fill(term, ident, fill, tag_handler, lat);
    }

    /** @brief Enum class used in the following class */
    enum class IonizedOrbital { Up, Down };

    /**
     * @brief Routine that builds the MPO representation of the destruction operator.
     * TODO: Check if it works also for the TwoU1 and SU2U1 symmetry.
     */
    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup>
    make_destroy_mpo(const Lattice& lat, const Model<Matrix, SymmGroup>& model, int iOrb, IonizedOrbital alphaOrBeta=IonizedOrbital::Up)
    {
      using tag_type = typename OPTable<Matrix, SymmGroup>::tag_type;
      std::vector<int> pos(1, iOrb);
      auto name = alphaOrBeta == IonizedOrbital::Up ? "destroy_up" : "destroy_down";
      tag_type opDestroy = model.get_operator_tag(name, lat.template get_prop<typename SymmGroup::subcharge>("type", pos[0]));
      tag_type ops_[1] = {opDestroy};
      std::vector<tag_type> ops(ops_, ops_+1);
      term_descriptor<typename Matrix::value_type> term = generate_mpo::arrange_operators(pos, ops, model.operators_table());
      std::vector<tag_type> identities, fillings;
      for (int iType = 0; iType < lat.getMaxType(); iType++) {
        identities.push_back(model.identity_matrix_tag(iType));
        fillings.push_back(model.filling_matrix_tag(iType));
      }
      return generate_mpo::make_1D_mpo(pos, ops, identities, fillings, model.operators_table(), lat);
    }

}

#endif
