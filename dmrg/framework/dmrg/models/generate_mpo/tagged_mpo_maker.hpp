/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef GENERATE_MPO_TAGGED_MPO_MAKER_H
#define GENERATE_MPO_TAGGED_MPO_MAKER_H

#include "dmrg/models/generate_mpo/utils.hpp"

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"

#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/chem_compression.h"

#include "dmrg/models/lattice.h"

#include <string>
#include <sstream>

namespace generate_mpo
{
    template<class Matrix, class SymmGroup>
    class TaggedMPOMaker
    {
        typedef typename Matrix::value_type scale_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
        typedef block_matrix<Matrix, SymmGroup> op_t;

        typedef Lattice::pos_t pos_t;
        typedef typename Operator_Tag_Term<Matrix, SymmGroup>::tag_type tag_type;
        typedef typename Operator_Tag_Term<Matrix, SymmGroup>::op_pair_t op_pair_t;
        typedef boost::tuple<std::size_t, std::size_t, tag_type, scale_type> tag_block;
        
    public:
        TaggedMPOMaker(pos_t length_, tag_type ident_
                      , boost::shared_ptr<OPTable<Matrix, SymmGroup> > tbl_)
        : length(length_)
        , used_dims(length)
        , ident(ident_)
        , tag_prempo(length)
        , maximum(2)
        , finalized(false)
        , leftmost_right(length)
        , operator_table(tbl_)
        {   
            for (size_t p = 0; p < length-1; ++p)
                tag_prempo[p].push_back(boost::make_tuple(std::size_t(0), std::size_t(0), ident, 1.));
        }
        
        void add_term(Operator_Tag_Term<Matrix, SymmGroup> const & term)
        {
            std::vector<op_pair_t> ops = term.operators;
            std::sort(ops.begin(), ops.end(), compare<std::pair<typename Lattice::pos_t, typename OPTable<Matrix, SymmGroup>::tag_type> >);
            
            vector<pos_t> positions;
            for (typename vector<op_pair_t>::const_iterator it = ops.begin(); it != ops.end(); ++it)
                positions.push_back(it->first);

            pos_t minp = *min_element(positions.begin(), positions.end());
            pos_t maxp = *max_element(positions.begin(), positions.end());
            
            std::size_t use_b = maximum++;
            
            vector<bool> done(length, false);
            for (typename vector<op_pair_t>::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
            {
                std::size_t first_use_b = (it->first == minp ? 0 : use_b);
                std::size_t second_use_b = (it->first == maxp ? 1 : use_b);
                assert(it->first < length);


                // apply scaling factor to first operator
                scale_type cur_scale = 1.;
                if (it == ops.begin())
                    cur_scale = term.scale;

                if (minp != maxp) { // bond term

                    /** tag_prempo filling **/
                    tag_prempo[it->first].push_back(boost::make_tuple(first_use_b, second_use_b, it->second, cur_scale));
                    used_dims[it->first].insert(use_b);

                } else { // site term

                    // retrieve the actual operator from the tag table
                    op_t current_op = (*operator_table)[it->second];
                    current_op *= term.scale;
                    site_terms[it->first] += current_op;
                }
                done[it->first] = true;
            }
            
            // put fill ops in between the nontrivial ops  (only needed for 1- and 2-terms)
            if (term.operators.size() <= 2) {
                op_t fill = (*operator_table)[term.fill_operator];
                for (pos_t p = minp; p <= maxp; ++p)
                    if (!done[p]) {
                        tag_prempo[p].push_back( boost::make_tuple(use_b, use_b, term.fill_operator, 1.) );
                        used_dims[p].insert(use_b);
                        done[p] = true;
                    }
            }
            
            leftmost_right = std::min(leftmost_right, maxp);
        }

        MPO<Matrix, SymmGroup> create_mpo()
        {
            if (!finalized) finalize(); 
            MPO<Matrix, SymmGroup> r(length);
            for (pos_t p = 1; p < length - 1; ++p)
                r[p] = as_bulk(tag_prempo[p]);
            r[0] = as_left(tag_prempo[0]);
            r[length-1] = as_right(tag_prempo[length-1]);
            
            return r;
        }

        /* Disabled for now
        MPO<Matrix, SymmGroup> create_compressed_mpo(Index<SymmGroup> const & phys, double cutoff)
        {
            if (!finalized) finalize();
            MPO<Matrix, SymmGroup> r = create_mpo();

            charge_sort(get_tag_prempo(), r);
            MPO<Matrix, SymmGroup> mpo_sorted = create_mpo();

            compressor<Matrix, SymmGroup> cpor(phys);
            MPO<Matrix, SymmGroup> mpo_out(length);
            cpor.compress(mpo_sorted, mpo_out, cutoff);

            return mpo_out;
        }
        */

        std::vector<std::vector<tag_block> > & get_tag_prempo()
        {
            return tag_prempo;
        }
        
    private:
        bool finalized;
        pos_t length;
        tag_type ident;
        vector<set<std::size_t> > used_dims;
        std::map<pos_t, op_t> site_terms;

        boost::shared_ptr<OPTable<Matrix, SymmGroup> > operator_table;
        vector<vector<tag_block> > tag_prempo;
        
        pos_t maximum, leftmost_right;

        void finalize()
        {
            for (typename std::map<pos_t, op_t>::const_iterator it = site_terms.begin();
                 it != site_terms.end(); ++it) {

                // TODO implement plus operation
                tag_type site_tag = operator_table->register_site_op(it->second);
                tag_prempo[it->first].push_back(boost::make_tuple(0, 1, site_tag, 1.));
            }

            for (pos_t p = leftmost_right + 1; p < length; ++p) {
                tag_prempo[p].push_back(boost::make_tuple(1, 1, ident, 1.));
            }

            for (typename vector<vector<tag_block> >::iterator it = tag_prempo.begin();
                 it + 1 != tag_prempo.end();
                 ++it)
                compress_on_bond(*it, *(it+1));

            finalized = true;
        }

        MPOTensor<Matrix, SymmGroup> as_bulk(vector<tag_block> const & ops)
        {
            pair<std::size_t, std::size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(rcd.first, rcd.second, ops, operator_table);
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_left(vector<tag_block> const & ops)
        {
            pair<std::size_t, std::size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(1, rcd.second, ops, operator_table);
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_right(vector<tag_block> const & ops)
        {
            pair<std::size_t, std::size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(rcd.first, 1, ops, operator_table);
            return r;
        }
    };

}

#endif
