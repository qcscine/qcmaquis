/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef GENERATE_MPO_TAGGED_MPO_MAKER_H
#define GENERATE_MPO_TAGGED_MPO_MAKER_H

#include "dmrg/models/generate_mpo/utils.hpp"

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"

#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mpo_ops.h"

#include "dmrg/models/lattice.h"
#include "dmrg/models/model.h"

#include <boost/bind.hpp>
#include <string>
#include <sstream>

namespace generate_mpo
{
    
    namespace detail {
        
        template <typename pos_t, typename tag_type, typename index_type>
        struct prempo_key {
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
    
    template <typename T, typename U>
    std::pair<T,U> to_pair(boost::tuple<T,U> const& t)
    {
        return std::make_pair( boost::get<0>(t), boost::get<1>(t) );
    }

    template<class Matrix, class SymmGroup>
    class TaggedMPOMaker
    {
        typedef typename Matrix::value_type scale_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
        typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;

        typedef Lattice::pos_t pos_t;
        typedef typename OperatorTagTerm<Matrix, SymmGroup>::tag_type tag_type;
        typedef typename OperatorTagTerm<Matrix, SymmGroup>::op_pair_t pos_op_type;
        typedef boost::tuple<std::size_t, std::size_t, tag_type, scale_type> tag_block;
        
        typedef ::term_descriptor<typename Matrix::value_type> term_descriptor;
        typedef std::vector<tag_type> tag_vec;
        
        typedef detail::prempo_key<pos_t, tag_type, index_type> prempo_key_type;
        typedef std::pair<tag_type, scale_type> prempo_value_type;
        // TODO: consider moving to hashmap
        typedef std::map<std::pair<prempo_key_type, prempo_key_type>, prempo_value_type> prempo_map_type;
        
        enum merge_kind {attach, detach};
        
    public:
        TaggedMPOMaker(Lattice const& lat_, Model<Matrix,SymmGroup> const& model)
        : lat(lat_)
        , length(lat.size())
        , tag_handler(model.operators_table())
        , prempo(length)
        , trivial_left(prempo_key_type::trivial_left)
        , trivial_right(prempo_key_type::trivial_right)
        , leftmost_right(length)
        , rightmost_left(0)
        , finalized(false)
        , verbose(true)
        , core_energy(0.)
        {
            for (size_t p = 0; p <= lat.maximum_vertex_type(); ++p)
            {
                identities.push_back(model.identity_matrix_tag(p));
                fillings.push_back(model.filling_matrix_tag(p));
                try { identities_full.push_back(model.get_operator_tag("ident_full", p)); }
                catch (std::runtime_error const & e) {}
            }

            typename Model<Matrix, SymmGroup>::terms_type const& terms = model.hamiltonian_terms();
            std::for_each(terms.begin(), terms.end(), boost::bind(&TaggedMPOMaker<Matrix,SymmGroup>::add_term, this, _1));
        }

        TaggedMPOMaker(Lattice const& lat_, tag_vec const & i_, tag_vec const & i_f_, tag_vec const & f_,
                       boost::shared_ptr<TagHandler<Matrix, SymmGroup> > th_, typename Model<Matrix, SymmGroup>::terms_type const& terms)
        : lat(lat_)
        , identities(i_)
        , identities_full(i_f_)
        , fillings(f_)
        , length(lat.size())
        , tag_handler(th_)
        , prempo(length)
        , trivial_left(prempo_key_type::trivial_left)
        , trivial_right(prempo_key_type::trivial_right)
        , leftmost_right(length)
        , rightmost_left(0)
        , finalized(false)
        , verbose(false)
        , core_energy(0.)
        {
            //for (size_t p = 0; p < length-1; ++p)
            //    prempo[p][make_pair(trivial_left,trivial_left)] = prempo_value_type(identities[lat.get_prop<int>("type",p)], 1.);
            
            std::for_each(terms.begin(), terms.end(), boost::bind(&TaggedMPOMaker<Matrix,SymmGroup>::add_term, this, _1));
        }
        
        void add_term(term_descriptor term)
        {
            std::sort(term.begin(), term.end(), pos_tag_lt());
            index_type nops = term.size();
            
            switch (nops) {
                case 1:
                    add_1term(term);
                    break;
                case 2:
                    add_2term(term);
                    break;
                case 3:
                    add_3term(term);
                    break;
                case 4:
                    add_4term(term);
                    break;
                default:
                    add_nterm(term); /// here filling has to be done manually
                    break;
            }
            
            leftmost_right = std::min(leftmost_right, boost::get<0>(*term.rbegin()));
            rightmost_left = std::max(rightmost_left, boost::get<0>(*term.begin()));
        }
                
        MPO<Matrix, SymmGroup> create_mpo()
        {
            if (!finalized) finalize();
            MPO<Matrix, SymmGroup> mpo; mpo.reserve(length);
            
            typedef std::map<prempo_key_type, index_type> index_map;
            typedef typename index_map::iterator index_iterator;
            index_map left;
            left[trivial_left] = 0;

            typedef SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type> spin_desc_t;
            std::vector<spin_desc_t> left_spins(1);
            
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
                    
                    index_type rr_dim = (p == length-1) ? 0 : rr->second;
                    pre_tensor.push_back( tag_block(ll->second, rr_dim, val.first, val.second) );
                }
                
                std::pair<index_type, index_type> rcd = rcdim(pre_tensor);

                std::vector<spin_desc_t> right_spins(rcd.second); 
                for (typename std::vector<tag_block>::const_iterator it = pre_tensor.begin(); it != pre_tensor.end(); ++it)
                {
                    spin_desc_t out_spin = couple(left_spins[boost::tuples::get<0>(*it)], tag_handler->get_op(boost::tuples::get<2>(*it)).spin());
                    index_type out_index = boost::tuples::get<1>(*it);
                    assert(right_spins[out_index].get() == 0 || right_spins[out_index].get() == out_spin.get());
                    right_spins[out_index] = out_spin;
                }

                if (p == 0)
                    mpo.push_back( MPOTensor<Matrix, SymmGroup>(1, rcd.second, pre_tensor, tag_handler->get_operator_table(), left_spins, right_spins) );
                else if (p == length - 1)
                    mpo.push_back( MPOTensor<Matrix, SymmGroup>(rcd.first, 1, pre_tensor, tag_handler->get_operator_table(), left_spins, right_spins) );
                else
                    mpo.push_back( MPOTensor<Matrix, SymmGroup>(rcd.first, rcd.second, pre_tensor, tag_handler->get_operator_table(), left_spins, right_spins) );
                if (verbose)
                    maquis::cout << "MPO Bond: " << rcd.second << std::endl;
                swap(left, right);
                swap(left_spins, right_spins);
            }
            
            mpo.setCoreEnergy(core_energy);
            return mpo;
        }
        
    private:
        void add_1term(term_descriptor const& term)
        {
            assert(term.size() == 1);
            
            /// Due to numerical instability: treat the core energy separately
            if (term.operator_tag(0) == identities[lat.get_prop<int>("type", term.position(0))])
                core_energy += term.coeff;

            else {
                /// retrieve the actual operator from the tag table
                op_t current_op = tag_handler->get_op(term.operator_tag(0));
                current_op *= term.coeff;
                site_terms[term.position(0)] += current_op;
            }
        }
        
        void add_2term(term_descriptor const& term)
        {
            assert(term.size() == 2);
            
            SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type > mpo_spin;
            prempo_key_type k1 = trivial_left;
            {
                int i = 0;
                mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
                prempo_key_type k2;
                k2.pos_op.push_back(to_pair(term[i+1]));
                k1 = insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), term.coeff), detach);
            }

            bool trivial_fill = !tag_handler->is_fermionic(term.operator_tag(1));
            // todo: check with long-range n_i*n_j                                  if spin > 0.5, need to use the full identity
            insert_filling(term.position(0)+1, term.position(1), k1, trivial_fill, (mpo_spin.get() > 1) ? term.full_identity : -1);
            {
                int i = 1;
                mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
                prempo_key_type k2 = trivial_right;
                insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), 1.), attach);
            }

            assert(mpo_spin.get() == 0); // H is a spin 0 operator
        }
        
        void add_3term(term_descriptor const& term)
        {
            assert(term.size() == 3);
            int nops = term.size();
            
            /// number of fermionic operators
            int nferm = 0;
            for (int i = 0; i < nops; ++i) {
                if (tag_handler->is_fermionic(term.operator_tag(i)))
                    nferm += 1;
            }

            SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type > mpo_spin;
            prempo_key_type k1 = trivial_left;
            std::vector<pos_op_type> ops_left;
            
            /// op_0
            {
                int i = 0;
                mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
                prempo_key_type k2;
                k2.pos_op.push_back(to_pair(term[i])); // k2: applied operator
                k1 = insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), 1.), attach);
                
                if (tag_handler->is_fermionic(term.operator_tag(i)))
                    nferm -= 1;
                bool trivial_fill = (nferm % 2 == 0);
                insert_filling(term.position(i)+1, term.position(i+1), k1, trivial_fill, (mpo_spin.get() > 1) ? term.full_identity : -1);
            }
            /// op_1
            {
                int i = 1;
                mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
                prempo_key_type k2;
                k2.pos_op.push_back(to_pair(term[i+1])); // k2: future operators
                k1 = insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), term.coeff), detach);
                
                if (tag_handler->is_fermionic(term.operator_tag(i)))
                    nferm -= 1;
                bool trivial_fill = (nferm % 2 == 0);
                insert_filling(term.position(i)+1, term.position(i+1), k1, trivial_fill, (mpo_spin.get() > 1) ? term.full_identity : -1);
            }
            /// op_2
            {
                int i = 2;
                mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
                insert_operator(term.position(i), make_pair(k1, trivial_right), prempo_value_type(term.operator_tag(i), 1.), attach);
            }

            assert(mpo_spin.get() == 0); // H is a spin 0 operator
        }
        
        void add_4term(term_descriptor const& term)
        {
            assert(term.size() == 4);
            int nops = term.size();
            
            /// number of fermionic operators
            int nferm = 0;
            for (int i = 0; i < nops; ++i) {
                if (tag_handler->is_fermionic(term.operator_tag(i)))
                    nferm += 1;
            }
            
            SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type > mpo_spin;
            prempo_key_type k1 = trivial_left;
            std::vector<pos_op_type> ops_left;
            
            /// op_0, op_1
            for (int i = 0; i < 2; ++i) {
                mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
                ops_left.push_back(to_pair(term[i])); prempo_key_type k2(ops_left);
                k1 = insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), 1.), attach);
                
                if (tag_handler->is_fermionic(term.operator_tag(i)))
                    nferm -= 1;
                bool trivial_fill = (nferm % 2 == 0);
                insert_filling(term.position(i)+1, term.position(i+1), k1, trivial_fill, (mpo_spin.get() > 1) ? term.full_identity : -1);
            }
            /// op_2
            {
                int i = 2;
                mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
                prempo_key_type k2;
                k2.pos_op.push_back(to_pair(term[3]));
                k1 = insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), term.coeff), detach);
                
                if (tag_handler->is_fermionic(term.operator_tag(i)))
                    nferm -= 1;
                bool trivial_fill = (nferm % 2 == 0);
                insert_filling(term.position(i)+1, term.position(i+1), k1, trivial_fill, (mpo_spin.get() > 1) ? term.full_identity : -1);
            }

            /// op_3
            {
                int i = 3;
                mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
                insert_operator(term.position(i), make_pair(k1, trivial_right), prempo_value_type(term.operator_tag(i), 1.), attach);
            }

            assert(mpo_spin.get() == 0); // H is a spin 0 operator
        }

        void add_nterm(term_descriptor const& term)
        {
            int nops = term.size();
            assert( nops > 2 );
            
            static index_type next_offset = 0;
            index_type current_offset = (next_offset++);
            
            prempo_key_type k1 = trivial_left;
            prempo_key_type k2(prempo_key_type::bulk_no_merge, current_offset);
            k2.pos_op.push_back( to_pair(term[nops-1]) );
            
            {
                int i = 0;
                insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), term.coeff), detach);
                k1 = k2;
                
                if (i < nops-1 && term.position(i)+1 != term.position(i+1))
                    throw std::runtime_error("for n > 4 operators filling is assumed to be done manually. the list of operators contains empty sites.");
            }

            
            for (int i = 1; i < nops; ++i) {
                if (i == nops-1)
                    k2 = trivial_right;
                
                insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), 1.), detach);
                
                if (i < nops-1 && term.position(i)+1 != term.position(i+1))
                    throw std::runtime_error("for n > 4 operators filling is assumed to be done manually. the list of operators contains empty sites.");
            }
            
        }

		void insert_filling(pos_t i, pos_t j, prempo_key_type k, bool trivial_fill, int custom_ident = -1)
		{
			for (; i < j; ++i) {
                tag_type use_ident = (custom_ident != -1) ? identities_full[lat.get_prop<int>("type",i)]
                                                          : identities[lat.get_prop<int>("type",i)];
                tag_type op = (trivial_fill) ? use_ident : fillings[lat.get_prop<int>("type",i)];
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
            if (merge_behavior == detach) {
                if (!match.second) {
                    std::pair<prempo_key_type, prempo_key_type> kk_max = kk;
                    kk_max.second.offset = std::numeric_limits<index_type>::max();

                    typename prempo_map_type::iterator highest_offset = prempo[p].upper_bound(kk_max);
                    highest_offset--;
                    kk.second.offset = highest_offset->first.second.offset + 1;
                    prempo[p].insert(highest_offset, make_pair(kk, val));
                }
            }
            else {
                // still slow, but seems never to be used
                while (!match.second && match.first->second != val) {
                    kk.second.offset += 1;
                    match = prempo[p].insert( make_pair(kk, val) );
                }
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

            // fill with ident from the begin
            for (size_t p = 0; p < rightmost_left; ++p)
                prempo[p][make_pair(trivial_left,trivial_left)] = prempo_value_type(identities[lat.get_prop<int>("type",p)], 1.);

            /// fill with ident until the end
            bool trivial_fill = true;
            insert_filling(leftmost_right+1, length, trivial_right, trivial_fill);

            finalized = true;
        }
        

    private:
        Lattice const& lat;

        tag_vec identities, identities_full, fillings;
        
        pos_t length;
        
        boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
        std::vector<prempo_map_type> prempo;
        prempo_key_type trivial_left, trivial_right;
        std::map<pos_t, op_t> site_terms;
        
        pos_t leftmost_right, rightmost_left;
        bool finalized, verbose;
        scale_type core_energy;
    };

}

#endif
