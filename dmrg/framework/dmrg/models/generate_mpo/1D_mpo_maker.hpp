/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2015-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef ONED_MPO_MAKER_H
#define ONED_MPO_MAKER_H

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
                      boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler)
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

            term.push_back( boost::make_tuple(pos_ops[opnr].first, product) );

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
                  boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler,
                  Lattice const & lat)
    {
        // after arrange operators, expand term to the full site-list

        typedef typename SymmGroup::subcharge sc;
        typedef Lattice::pos_t pos_t;
        typedef term_descriptor<typename Matrix::value_type> term_descriptor;
        typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
        typename MPOTensor<Matrix, SymmGroup>::prempo_t prempo;

        using boost::tuples::get;

        MPO<Matrix, SymmGroup> ret(lat.size());

        bool carry_sign = false;
        for (pos_t p = 0; p < lat.size(); ++p)
        {
            // check if there is a non-trivial operator on site p
            typename term_descriptor::const_iterator
                it = std::find_if(term.begin(), term.end(),
                    boost::bind(detail::get<0, typename term_descriptor::value_type>, boost::lambda::_1) == p);

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
                boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler,
                Lattice const & lat)
    {
        term_descriptor<typename Matrix::value_type> term = arrange_operators(positions, operators, tag_handler);
        return sign_and_fill(term, ident, fill, tag_handler, lat);
    }

}

#endif
