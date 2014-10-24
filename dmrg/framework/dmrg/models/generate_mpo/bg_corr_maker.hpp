/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef BG_CORR_MAKER_H
#define BG_CORR_MAKER_H

#include "dmrg/models/generate_mpo/corr_maker.hpp"

namespace generate_mpo
{
    
    template <class Matrix, class SymmGroup>
    std::vector<std::pair<Lattice::pos_t, std::pair<std::vector<block_matrix<Matrix, SymmGroup> >, bool> > >
    arrange_operators(std::vector<Lattice::pos_t> const & positions,
                      std::vector<std::pair<std::vector<block_matrix<Matrix, SymmGroup> >, bool> > const & operators)

    {
        // input: list of positions and operators
        // output: list of (position, operator)-pairs, sorted, unique positions with operators multiplied
        typedef Lattice::pos_t pos_t;
        typedef block_matrix<Matrix, SymmGroup> op_t;
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

    template<class Matrix, class SymmGroup>
    class BgCorrMaker : public CorrMakerBase<Matrix, SymmGroup>
    {
        typedef CorrMakerBase<Matrix, SymmGroup> base;
        typedef tag_detail::tag_type tag_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
        typedef typename base::block block;
        typedef boost::tuple<size_t, size_t, string> tag;

        typedef Lattice::pos_t pos_t;
        typedef block_matrix<Matrix, SymmGroup> op_t;
        typedef std::pair<std::vector<op_t>, bool> site_ops_t;
        typedef std::pair<pos_t, site_ops_t> pos_op_t;
        
    public:
        BgCorrMaker(Lattice const& lat_,
                  const std::vector<op_t> & ident_,
                  const std::vector<op_t> & fill_,
                  std::vector<site_ops_t> const & ops,
                  std::vector<pos_t> const & ref,
                  bool idiag = false)
        : lat(lat_)
        , prempo(lat.size())
        , tags(lat.size())
        , identities(ident_.size())
        , fillings(fill_.size())
        , op_tags(ops.size())
        , incl_diag(idiag)
        {
            /// register operators
            for (int type=0; type<ident_.size(); ++type)
                identities[type] = tag_handler.register_op(ident_[type], tag_detail::bosonic);
            for (int type=0; type<fill_.size(); ++type)
                fillings[type] = tag_handler.register_op(fill_[type], tag_detail::bosonic);
            
            for (size_t n=0; n<ops.size(); ++n) {
                op_tags[n].resize(ops[n].first.size());
                std::vector<tag_type> & tops = op_tags[n];
                for (int type=0; type<tops.size(); ++type)
                    tops[type] = tag_handler.register_op(ops[n].first[type],
                                    ops[n].second ? tag_detail::fermionic : tag_detail::bosonic);
            }
            if(lat.size() < 2)
                throw std::runtime_error("Need lattice size > 1 for correlation measurements\n");

            if(ref.size() != ops.size()-1)
                throw std::runtime_error("CorrMaker was called with wrong number of fixed operator positions\n");

            std::copy(ref.begin(), ref.end()-1, std::back_inserter(background_pos));
            phase = set_base_phase(ref);
            labels.resize(lat.size() - *ref.rbegin() - 1 + (int)incl_diag);

            // Handle operators on identical positions
            std::vector<pos_op_t> pos_ops = arrange_operators(ref, ops);

            for (int i=0; i < ref.size()-1; ++i)
                for (int j=0; j < pos_ops.size(); ++j)
                {
                    if (ref[i] == pos_ops[j].first)
                        for (size_t type=0; type < ops[0].first.size(); ++type)
                            op_tags[i][type] = tag_handler.register_op(pos_ops[j].second.first[type],
                                                                       pos_ops[j].second.second ? tag_detail::fermionic : tag_detail::bosonic);
                    break;
                }

            make_prempo(*ref.rbegin());
        }
        
        MPO<Matrix, SymmGroup> create_mpo()
        {
            boost::shared_ptr<OPTable<Matrix, SymmGroup> > tbl = tag_handler.get_operator_table();
            MPO<Matrix, SymmGroup> r(prempo.size());
            for (size_t p = 1; p < prempo.size(); ++p)
                r[p] = base::as_bulk(prempo[p], tbl);
            r[0] = base::as_left(prempo[0], tbl);
            
            return r;
        }
        
        std::string description () const
        {
            std::ostringstream ss;
        	for (size_t p = 0; p < prempo.size(); ++p)
            {
                ss << "Site: " << p << std::endl;
                for (typename vector<tag>::const_iterator it = tags[p].begin(); it != tags[p].end(); ++it)
                    ss << "    " << get<0>(*it) << " " << get<1>(*it) << " " << get<2>(*it) << std::endl;
            }
        	return ss.str();
        }
        
        vector<vector<pos_t> > const& numeric_labels() { return labels; }
        
    private:
        Lattice const& lat;
        TagHandler<Matrix, SymmGroup> tag_handler;

        vector<vector<block> > prempo;
        vector<vector<tag> > tags;
        vector<vector<pos_t> > labels;
        
        std::vector<tag_type> identities, fillings;
        vector<std::vector<tag_type> > op_tags;
        std::vector<pos_t> background_pos;
        bool incl_diag;
        int phase; // needed for background, if present

        void make_prempo(pos_t start)
        {
            bool start_sign = false;
            for(pos_t p = 0; p < start; ++p)
                start_sign = insert_filling(0, 0, p, start_sign);

            bool current_sign = insert_operator(0, 0, start, (*(op_tags.rbegin()+1))[lat.get_prop<int>("type", start)], start_sign);

            for(pos_t branch = start+1; branch < lat.size(); ++branch)
            {
                size_t branch_index = lat.size() - branch - 1;
                if (branch > start+1)
                    current_sign = insert_filling(0, 0, branch-1, current_sign);
                bool branch_sign = insert_operator(0, branch_index, branch,
                                                   (*op_tags.rbegin())[lat.get_prop<int>("type", branch)], current_sign);
                for(pos_t p2 = branch+1; p2 < lat.size(); ++p2)
                    branch_sign = insert_filling(branch_index, branch_index, p2, branch_sign);

                labels[branch_index] = background_pos;
                labels[branch_index].push_back(start);
                labels[branch_index].push_back(branch);
            }

            if(incl_diag)
            {
                tag_type op_tag;
                typename Matrix::value_type scale;
                size_t branch_index = lat.size() - start - 1;
			    boost::tie(op_tag, scale) = tag_handler.get_product_tag((*op_tags.rbegin())[lat.get_prop<int>("type", start)],
                                                                        (*(op_tags.rbegin()+1))[lat.get_prop<int>("type", start)]);

                bool branch_sign = insert_operator(0, branch_index, start, op_tag, start_sign, scale);
                for(pos_t p2 = start+1; p2 < lat.size(); ++p2)
                    branch_sign = insert_filling(branch_index, branch_index, p2, branch_sign);

                labels[branch_index] = background_pos;
                labels[branch_index].push_back(start);
                labels[branch_index].push_back(start);
            }
        }

        bool insert_filling(index_type b1, index_type b2, pos_t p, bool sign_in)
        {
            // check if a background operator is present
            std::vector<pos_t>::const_iterator it = std::find(background_pos.begin(), background_pos.end(), p);
            if( it != background_pos.end() ) {
                int bg_index = it - background_pos.begin();
                return insert_operator_(b1, b2, p, op_tags[bg_index][lat.get_prop<int>("type", p)], sign_in);
            }

            if (sign_in)
            prempo[p].push_back(boost::make_tuple(b1, b2, fillings[lat.get_prop<int>("type", p)], 1.));
            else
            prempo[p].push_back(boost::make_tuple(b1, b2, identities[lat.get_prop<int>("type", p)], 1.));
            return sign_in;
        }

        bool insert_operator(index_type b1, index_type b2, pos_t p,
                             tag_type opt, bool sign_in, typename Matrix::value_type scale_in=1.)
        {
            // check if a background operator is present
            std::vector<pos_t>::const_iterator it = std::find(background_pos.begin(), background_pos.end(), p);
            if( it != background_pos.end() ) {
                int bg_index = it - background_pos.begin();
                typename Matrix::value_type scale;
                tag_type op_tag;

                // multiply op with background operator
			    boost::tie(op_tag, scale) = tag_handler.get_product_tag(opt, op_tags[bg_index][lat.get_prop<int>("type", p)]);

                return insert_operator_(b1, b2, p, op_tag, sign_in, scale * scale_in);
            }

            return insert_operator_(b1, b2, p, opt, sign_in, scale_in);
        }

        bool insert_operator_(index_type b1, index_type b2, pos_t p, tag_type op,
                              bool sign_in, typename Matrix::value_type scale_in=1.)
        {
            // Here, a phase factor stemming from the background operators (if present)
            // is applied. Condition below = indirect way to identify the last nontriv op per branch
            if( b1 != b2 || (p == lat.size()-1 && b1==0 && b2==0) )
                scale_in *= final_phase(p);

            bool sign_out = (sign_in != tag_handler.is_fermionic(op));
            typename Matrix::value_type scale_loc = 1.;
            tag_type fill = fillings[lat.get_prop<int>("type", p)];
            if (sign_out)
			    boost::tie(op, scale_loc) = tag_handler.get_product_tag(fill, op);

            prempo[p].push_back(boost::make_tuple(b1, b2, op, scale_loc * scale_in));
            return sign_out;
        }

        int set_base_phase(std::vector<pos_t> const & pos)
        {
            int inv_count=0, n=pos.size();
            for(pos_t c1 = 0; c1 < n - 1; c1++)
                for(pos_t c2 = c1+1; c2 < n; c2++)
                    if(pos[c1] > pos[c2]) inv_count++; 
            return inv_count;
        }

        int final_phase(pos_t p)
        {
            int cnt=0;
            for(size_t k = 0; k < background_pos.size(); ++k)
                if (p < background_pos[k]) cnt++;

            return 1 - 2 * ((cnt + phase) % 2);
        }
    };
}

#endif
