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
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/mp_tensors/chem_compression.h"

#include "dmrg/models/lattice.h"

#include <string>
#include <sstream>

namespace generate_mpo
{
    
    namespace detail {
        
        template <typename pos_t, typename tag_type, typename index_type>
        struct prempo_key {
//            typedef Lattice::pos_t pos_t;
//            typedef tag_detail::tag_type tag_type;
//            typedef unsigned index_type;
            typedef std::pair<pos_t, tag_type> pos_op_type;
            enum kind_type {trivial_left, bulk, bulk_no_merge, trivial_right};
            
            kind_type kind;
            std::vector<pos_op_type> pos_op;
            index_type offset;
            
            prempo_key(kind_type k_=bulk, index_type o_=0) : kind(k_), offset(o_) { }
            prempo_key(std::vector<pos_op_type> const& po_, index_type o_=0) : kind(bulk), pos_op(po_), offset(o_) { }
            
            bool operator==(prempo_key const& lhs) const
            {
                if (kind != lhs.kind)
                    return false;
                if (kind == trivial_left)
                    return true;
                if (kind == trivial_right)
                    return true;
                
                return (pos_op == lhs.pos_op) && (offset == lhs.offset);
            }
            
            bool operator<(prempo_key const& lhs) const
            {
                if (kind != lhs.kind) return kind < lhs.kind;
                return (pos_op == lhs.pos_op) ? offset < lhs.offset : pos_op < lhs.pos_op;
            }
        };
        
        
    }
    
    template<class Matrix, class SymmGroup>
    class TaggedMPOMaker
    {
        typedef typename Matrix::value_type scale_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
        typedef block_matrix<Matrix, SymmGroup> op_t;

        typedef Lattice::pos_t pos_t;
        typedef typename Operator_Tag_Term<Matrix, SymmGroup>::tag_type tag_type;
        typedef typename Operator_Tag_Term<Matrix, SymmGroup>::op_pair_t pos_op_type;
        typedef boost::tuple<std::size_t, std::size_t, tag_type, scale_type> tag_block;
        
        typedef detail::prempo_key<pos_t, tag_type, index_type> prempo_key_type;
        typedef std::pair<tag_type, scale_type> prempo_value_type;
        // TODO: consider moving to hashmap
        typedef std::map<std::pair<prempo_key_type, prempo_key_type>, prempo_value_type> prempo_map_type;
        
        enum merge_kind {attach, detach};
        
    public:
        TaggedMPOMaker(pos_t length_, tag_type ident_
                      , boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tbl_)
        : length(length_)
        , ident(ident_)
        , prempo(length)
        , trivial_left(prempo_key_type::trivial_left)
        , trivial_right(prempo_key_type::trivial_right)
        , leftmost_right(length)
        , tag_handler(tbl_)
        , finalized(false)
        {
            for (size_t p = 0; p < length-1; ++p)
                prempo[p][make_pair(trivial_left,trivial_left)] = prempo_value_type(ident, 1.);
        }
        
        void add_term(Operator_Tag_Term<Matrix, SymmGroup> const & term)
        {
            std::vector<pos_op_type> ops = term.operators;
            std::sort(ops.begin(), ops.end(), compare<pos_op_type>);
            index_type nops = ops.size();
            
            switch (nops) {
                case 1:
                    add_1term(ops, term.scale, term.fill_operator);
                    break;
                case 2:
                    add_2term(ops, term.scale, term.fill_operator);
                    break;
                case 3:
                    add_3term(ops, term.scale, term.fill_operator);
                    break;
                case 4:
                    add_4term(ops, term.scale, term.fill_operator);
                    break;
                default:
                    add_nterm(ops, term.scale); /// here filling has to be done manually
                    break;
            }
            
            leftmost_right = std::min(leftmost_right, ops.rbegin()->first);
        }
                
        MPO<Matrix, SymmGroup> create_mpo()
        {
            if (!finalized) finalize();
            MPO<Matrix, SymmGroup> mpo; mpo.reserve(length);
            
            typedef std::map<prempo_key_type, index_type> index_map;
            typedef typename index_map::iterator index_iterator;
            index_map left;
            left[trivial_left] = 0;
            
            for (pos_t p = 0; p < length; ++p) {
                std::vector<tag_block> pre_tensor; pre_tensor.reserve(prempo[p].size());
                
                index_map right;
                index_type r = 2;
                for (typename prempo_map_type::const_iterator it = prempo[p].begin();
                     it != prempo[p].end(); ++it) {
                    prempo_key_type const& k1 = it->first.first;
                    prempo_key_type const& k2 = it->first.second;
                    prempo_value_type const& val = it->second;
                    
                    index_iterator ll = left.find(k1);
                    if (ll == left.end())
                        throw std::runtime_error("k1 not found!");
                    
                    index_iterator rr = right.find(k2);
                    if (k2 == trivial_left && rr == right.end())
                        boost::tie(rr, boost::tuples::ignore) = right.insert( make_pair(k2, 0) );
                    else if (k2 == trivial_right && rr == right.end())
                        boost::tie(rr, boost::tuples::ignore) = right.insert( make_pair(k2, 1) );
                    else if (rr == right.end())
                        boost::tie(rr, boost::tuples::ignore) = right.insert( make_pair(k2, r++) );
                    
                    pre_tensor.push_back( tag_block(ll->second, rr->second, val.first, val.second) );
                }
                
                std::pair<std::size_t, std::size_t> rcd = rcdim(pre_tensor);
                if (p == 0)
                    mpo.push_back( MPOTensor<Matrix, SymmGroup>(1, rcd.second, pre_tensor, tag_handler->get_operator_table()) );
                else if (p == length - 1)
                    mpo.push_back( MPOTensor<Matrix, SymmGroup>(rcd.first, 1, pre_tensor, tag_handler->get_operator_table()) );
                else
                    mpo.push_back( MPOTensor<Matrix, SymmGroup>(rcd.first, rcd.second, pre_tensor, tag_handler->get_operator_table()) );

                swap(left, right);
            }
            
            follow_and_print_terms(mpo, -1, 0, 0);
            return mpo;
        }
        
    private:
        void add_1term(std::vector<pos_op_type> const& ops, scale_type matrix_element, tag_type fill)
        {
            assert(ops.size() == 1);
            
            /// retrieve the actual operator from the tag table
            // TODO implement plus operation
            op_t current_op = tag_handler->get_op(ops[0].second);
            current_op *= matrix_element;
            site_terms[ops[0].first] += current_op;
        }
        
        void add_2term(std::vector<pos_op_type> const& ops, scale_type matrix_element, tag_type fill)
        {
            assert(ops.size() == 2);
            
            prempo_key_type k1 = trivial_left;
            {
                int i = 0;
                prempo_key_type k2;
                k2.pos_op.push_back(ops[i+1]);
                k1 = insert_operator(ops[i].first, make_pair(k1, k2), prempo_value_type(ops[i].second, matrix_element), detach);
            }
            insert_filling(ops[0].first+1, ops[1].first, k1, fill);
            {
                int i = 1;
                prempo_key_type k2 = trivial_right;
                insert_operator(ops[i].first, make_pair(k1, k2), prempo_value_type(ops[i].second, 1.), detach);
            }
}
        
        void add_3term(std::vector<pos_op_type> const& ops, scale_type matrix_element, tag_type fill)
        {
            assert(ops.size() == 3);
            int nops = ops.size();
            
            /// number of fermionic operators
            int nferm = 0;
            for (int i = 0; i < nops; ++i) {
                if (tag_handler->is_fermionic(ops[i].second))
                    nferm += 1;
            }
            
            prempo_key_type k1 = trivial_left;
            std::vector<pos_op_type> ops_left;
            
            /// op_0
            {
                int i = 0;
                prempo_key_type k2;
                k2.pos_op.push_back(ops[i]); // k2: applied operator
                k1 = insert_operator(ops[i].first, make_pair(k1, k2), prempo_value_type(ops[i].second, 1.), attach);
                
                if (tag_handler->is_fermionic(ops[i].second))
                    nferm -= 1;
                tag_type curr_fill = (nferm % 2 == 0) ? ident : fill;
                insert_filling(ops[i].first+1, ops[i+1].first, k1, curr_fill);
            }
            /// op_1
            {
                int i = 1;
                prempo_key_type k2;
                k2.pos_op.push_back(ops[i+1]); // k2: future operators
                k1 = insert_operator(ops[i].first, make_pair(k1, k2), prempo_value_type(ops[i].second, matrix_element), detach);
                
                if (tag_handler->is_fermionic(ops[i].second))
                    nferm -= 1;
                tag_type curr_fill = (nferm % 2 == 0) ? ident : fill;
                insert_filling(ops[i].first+1, ops[i+1].first, k1, curr_fill);
            }
            /// op_2
            {
                int i = 2;
                insert_operator(ops[i].first, make_pair(k1, trivial_right), prempo_value_type(ops[i].second, 1.), detach);
            }
        }
        
        void add_4term(std::vector<pos_op_type> const& ops, scale_type matrix_element, tag_type fill)
        {
            assert(ops.size() == 4);
            int nops = ops.size();
            
            /// number of fermionic operators
            int nferm = 0;
            for (int i = 0; i < nops; ++i) {
                if (tag_handler->is_fermionic(ops[i].second))
                    nferm += 1;
            }
            
            prempo_key_type k1 = trivial_left;
            std::vector<pos_op_type> ops_left;
            
            /// op_0, op_1
            for (int i = 0; i < 2; ++i) {
                ops_left.push_back(ops[i]); prempo_key_type k2(ops_left);
                k1 = insert_operator(ops[i].first, make_pair(k1, k2), prempo_value_type(ops[i].second, 1.), attach);
                
                if (tag_handler->is_fermionic(ops[i].second))
                    nferm -= 1;
                tag_type curr_fill = (nferm % 2 == 0) ? ident : fill;
                insert_filling(ops[i].first+1, ops[i+1].first, k1, curr_fill);
            }
            /// op_2
            {
                int i = 2;
                prempo_key_type k2;
                k2.pos_op.push_back(ops[3]);
                k1 = insert_operator(ops[i].first, make_pair(k1, k2), prempo_value_type(ops[i].second, matrix_element), detach);
                
                if (tag_handler->is_fermionic(ops[i].second))
                    nferm -= 1;
                tag_type curr_fill = (nferm % 2 == 0) ? ident : fill;
                insert_filling(ops[i].first+1, ops[i+1].first, k1, curr_fill);
            }

            /// op_3
            {
                int i = 3;
                insert_operator(ops[i].first, make_pair(k1, trivial_right), prempo_value_type(ops[i].second, 1.), detach);
            }
        }

        void add_nterm(std::vector<pos_op_type> const& ops, scale_type matrix_element)
        {
            int nops = ops.size();
            assert( nops > 2 );
            
            static index_type next_offset = 0;
            index_type current_offset = (next_offset++);
            
            prempo_key_type k1 = trivial_left;
            prempo_key_type k2(prempo_key_type::bulk_no_merge, current_offset);
            k2.pos_op.push_back( ops[nops-1] );
            
            {
                int i = 0;
                insert_operator(ops[i].first, make_pair(k1, k2), prempo_value_type(ops[i].second, matrix_element), detach);
                k1 = k2;
                
                if (i < nops-1 && ops[i].first+1 != ops[i+1].first)
                    throw std::runtime_error("for n > 4 operators filling is assumed to be done manually. the list of operators contains empty sites.");
            }

            
            for (int i = 1; i < nops; ++i) {
                if (i == nops-1)
                    k2 = trivial_right;
                
                insert_operator(ops[i].first, make_pair(k1, k2), prempo_value_type(ops[i].second, 1.), detach);
                
                if (i < nops-1 && ops[i].first+1 != ops[i+1].first)
                    throw std::runtime_error("for n > 4 operators filling is assumed to be done manually. the list of operators contains empty sites.");
            }
            
        }

		void insert_filling(pos_t i, pos_t j, prempo_key_type k, tag_type op)
		{
			for (; i < j; ++i) {
				std::pair<typename prempo_map_type::iterator,bool> ret = prempo[i].insert( make_pair(make_pair(k,k), prempo_value_type(op, 1.)) );
				if (!ret.second && ret.first->second.first != op)
					throw std::runtime_error("Pre-existing term at site "+boost::lexical_cast<std::string>(i)
					                         + ". Needed "+boost::lexical_cast<std::string>(op)
					                         + ", found "+boost::lexical_cast<std::string>(ret.first->second.first));
			}
		}

		prempo_key_type insert_operator(pos_t p, std::pair<prempo_key_type, prempo_key_type> kk, prempo_value_type val,
                                        merge_kind merge_behavior=detach)
		{
			/// merge_behavior == detach: a new branch will be created, in case op already exist, an offset is used
			/// merge_behavior == attach: if operator tags match, keep the same branch
            std::pair<typename prempo_map_type::iterator,bool> match = prempo[p].insert( make_pair(kk, val) );
			while (!match.second && (merge_behavior==detach || match.first->second != val)) {
				kk.second.offset += 1;
				match = prempo[p].insert( make_pair(kk, val) );
			}
			return kk.second;
		}
		
        void finalize()
        {
            /// site terms
            std::pair<prempo_key_type,prempo_key_type> kk = make_pair(trivial_left,trivial_right);
            for (typename std::map<pos_t, op_t>::const_iterator it = site_terms.begin();
                 it != site_terms.end(); ++it) {
                tag_type site_tag = tag_handler->register_op(it->second, tag_detail::bosonic);
				std::pair<typename prempo_map_type::iterator,bool> ret;
                ret = prempo[it->first].insert( make_pair( kk, prempo_value_type(site_tag,1.) ) );
                if (!ret.second)
                    throw std::runtime_error("another site term already existing!");
            }

            /// fill with ident until the end
            insert_filling(leftmost_right+1, length, trivial_right, ident);

            finalized = true;
        }
        

    private:
        pos_t length;
        tag_type ident;
        
        boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
        std::vector<prempo_map_type> prempo;
        prempo_key_type trivial_left, trivial_right;
        std::map<pos_t, op_t> site_terms;
        
        pos_t leftmost_right;
        bool finalized;
    };

}

#endif
