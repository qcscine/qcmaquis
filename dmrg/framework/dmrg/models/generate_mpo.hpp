/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef GENERATE_MPO_H
#define GENERATE_MPO_H

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"

#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/chem_compression.h"

#include "dmrg/models/lattice.h"

#include <string>
#include <sstream>

#include <boost/bind.hpp>

namespace generate_mpo
{
	template<class Matrix, class SymmGroup>
	struct Operator_Tag_Term
	{
		typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
        typedef typename Lattice::pos_t pos_t;
		typedef std::pair<pos_t, tag_type> op_pair_t;
        
		std::vector<op_pair_t> operators;
		tag_type fill_operator;
        typename Matrix::value_type scale;
        bool with_sign;
        
        Operator_Tag_Term() : scale(1.), with_sign(false) {}
	};
    
	template<class Matrix, class SymmGroup>
	struct Operator_Term
	{
		typedef block_matrix<Matrix, SymmGroup> op_t;
        typedef typename Lattice::pos_t pos_t;
		typedef std::pair<pos_t, op_t> op_pair_t;
        
		std::vector<op_pair_t> operators;
		op_t fill_operator;
        bool with_sign;
        
        Operator_Term() : with_sign(false) {}
        
        void canonical_order() // TODO: check and fix for fermions
        {
            std::sort(operators.begin(), operators.end(),
                      boost::bind(&op_pair_t::first, _1) <
                      boost::bind(&op_pair_t::first, _2));
        }
        
        bool operator< (Operator_Term const & rhs) const
        {
            if (operators[0].first == rhs.operators[0].first)
                return operators.size() >= rhs.operators.size();
            return operators[0].first < rhs.operators[0].first;
        }

        bool site_match (Operator_Term const & rhs) const
        {
            if (operators.size() == rhs.operators.size())
            {
                bool ret = true;
                for (std::size_t p=0; p<operators.size() && ret; ++p)
                    ret = (operators[p].first == rhs.operators[p].first);
                return ret;
            } else if (operators.size() == 2 && rhs.operators.size() == 1)
                return (operators[0].first == rhs.operators[0].first || operators[1].first == rhs.operators[0].first);
            else if (operators.size() == 1 && rhs.operators.size() == 2)
                return (operators[0].first == rhs.operators[0].first || operators[0].first == rhs.operators[1].first);
            else
            {
                throw std::runtime_error("site_match not implemented for this type of operator." );
                return false;
            }
                
        }
        
        bool overlap (Operator_Term const & rhs) const
        {
        	return !( (operators.rbegin()->first < rhs.operators.begin()->first) || (rhs.operators.rbegin()->first < operators.begin()->first) );
        }

	};
   
    using namespace std;
    using namespace boost::tuples;

    inline size_t next_free(vector<size_t> const & out_taken,
                            vector<size_t> const & in_taken)
    {
        for (size_t k = 0; true; ++k)
        {
            if (count(out_taken.begin(), out_taken.end(), k) == 0 &&
                count(in_taken.begin(), in_taken.end(), k) == 0)
                return k;
        }
    }
    
    inline size_t next_free(set<size_t> const & s)
    {
        for (size_t k = 2; true; ++k)
            if (s.count(k) == 0)
                return k;
    }
    
    template<class Vector>
    void compress_on_bond(Vector & pm1, Vector & pm2)
    {
        std::set<size_t> bond_used_dims;
        for (typename Vector::iterator it = pm1.begin(); it != pm1.end(); ++it)
            if (get<1>(*it) > 1)
                bond_used_dims.insert(get<1>(*it));
        for (typename Vector::iterator it = pm2.begin(); it != pm2.end(); ++it)
            if (get<0>(*it) > 1)
                bond_used_dims.insert(get<0>(*it));
        
        std::map<size_t, size_t> compression_map;
        size_t c = 2;
        for (set<size_t>::iterator it = bond_used_dims.begin();
             it != bond_used_dims.end(); ++it)
            compression_map[*it] = c++;
        
        for (typename Vector::iterator it = pm1.begin(); it != pm1.end(); ++it)
            if (compression_map.count(get<1>(*it)) > 0)
                get<1>(*it) = compression_map[get<1>(*it)];
        for (typename Vector::iterator it = pm2.begin(); it != pm2.end(); ++it)
            if (compression_map.count(get<0>(*it)) > 0)
                get<0>(*it) = compression_map[get<0>(*it)];
    }

    template<class Vector>
    std::pair<size_t, size_t> rcdim(Vector const & pm)
    {
        std::list<size_t> l, r;
        for (typename Vector::const_iterator it = pm.begin(); it != pm.end(); ++it) {
            l.push_back( get<0>(*it) );
            r.push_back( get<1>(*it) );
        }
        
        return make_pair(*max_element(l.begin(), l.end())+1,
                         *max_element(r.begin(), r.end())+1);
    }
    
    template<class Pair>
    bool compare(Pair const & p1, Pair const & p2)
    {
        return p1.first < p2.first;
    }

    template<class Matrix, class SymmGroup>
    class MPOMaker
    {
        typedef block_matrix<Matrix, SymmGroup> op_t;
        typedef boost::tuple<size_t, size_t, op_t> block;
        typedef vector<
        pair<
        block_matrix<Matrix, SymmGroup>,
        block_matrix<Matrix, SymmGroup>
        >
        > op_pairs;
        
    public:
        MPOMaker(std::size_t length_, const block_matrix<Matrix, SymmGroup> & ident_)
        : length(length_)
        , used_dims(length)
        , ident(ident_)
        , prempo(length)
        , maximum(2)
        , finalized(false)
        , leftmost_right(length)
        {   
            for (size_t p = 0; p < length; ++p)
            {
                if (p+1 < length)
                    prempo[p].push_back(boost::make_tuple(std::size_t(0), std::size_t(0), ident));
            }
        }
        
        void add_term(Operator_Term<Matrix, SymmGroup> const & term)
        {
            // TODO: removed const&, because of sorting (non-const operation)
            std::vector<std::pair<typename Lattice::pos_t, op_t> > ops = term.operators;
            
            std::sort(ops.begin(), ops.end(), compare<std::pair<typename Lattice::pos_t, op_t> >);
            
            vector<size_t> positions;
            for (typename vector<pair<typename Lattice::pos_t, op_t> >::const_iterator
                 it = ops.begin();
                 it != ops.end(); ++it)
                positions.push_back( it->first );
            size_t minp = *min_element(positions.begin(), positions.end());
            size_t maxp = *max_element(positions.begin(), positions.end());
            
            size_t use_b = maximum++;
            
            vector<bool> done(length, false);
            for (typename vector<pair<typename Lattice::pos_t, op_t> >::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
            {
                size_t first_use_b = (it->first == minp ? 0 : use_b);
                size_t second_use_b = (it->first == maxp ? 1 : use_b);
                assert( it->first < prempo.size() );
                if (minp != maxp) { // bond term
                    prempo[it->first].push_back(boost::make_tuple(first_use_b, second_use_b, it->second));
                    used_dims[it->first].insert(use_b);
                } else // site term
                    site_terms[it->first] += it->second;
                done[it->first] = true;
            }
            
            for (size_t p = minp; p <= maxp; ++p)
                if (!done[p]) {
                    prempo[p].push_back( boost::make_tuple(use_b, use_b, term.fill_operator));
                    used_dims[p].insert(use_b);
                    done[p] = true;
                }
            
            leftmost_right = std::min(leftmost_right, maxp);
        }
        
        MPO<Matrix, SymmGroup> create_mpo()
        {
            if (!finalized) finalize(); 
            MPO<Matrix, SymmGroup> r(length);
            for (size_t p = 1; p < length - 1; ++p)
                r[p] = as_bulk(prempo[p]);
            r[0] = as_left(prempo[0]);
            r[length-1] = as_right(prempo[length-1]);
            
            return r;
        }

        MPO<Matrix, SymmGroup> create_compressed_mpo(Index<SymmGroup> const & phys, double cutoff)
        {
            if (!finalized) finalize();
            MPO<Matrix, SymmGroup> r(length);
            for (size_t p = 1; p < length - 1; ++p)
                r[p] = as_bulk(prempo[p]);
            r[0] = as_left(prempo[0]);
            r[length-1] = as_right(prempo[length-1]);

            charge_sort(get_prempo(), r);
            MPO<Matrix, SymmGroup> mpo_sorted = create_mpo();

            compressor<Matrix, SymmGroup> cpor(phys);
            MPO<Matrix, SymmGroup> mpo_out(length);
            cpor.compress(mpo_sorted, mpo_out, cutoff);

            return mpo_out;
        }

        std::vector<std::vector<block> > & get_prempo()
        {
            return prempo;
        }
        
    private:
        bool finalized;
        std::size_t length;
        block_matrix<Matrix, SymmGroup> ident;
        vector<set<size_t> > used_dims;
        vector<vector<block> > prempo;
        std::map<std::size_t, op_t> site_terms;
        
        size_t maximum, leftmost_right;

        void finalize()
        {
            for (typename std::map<std::size_t, op_t>::const_iterator it = site_terms.begin();
                 it != site_terms.end(); ++it)
                prempo[it->first].push_back( boost::make_tuple(0, 1, it->second) );

            for (size_t p = leftmost_right + 1; p < length; ++p)
                prempo[p].push_back( boost::make_tuple(1, 1, ident) );

            for (typename vector<vector<block> >::iterator it = prempo.begin();
                 it + 1 != prempo.end();
                 ++it)
                compress_on_bond(*it, *(it+1));

            finalized = true;
        }
        
        MPOTensor<Matrix, SymmGroup> as_bulk(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(rcd.first, rcd.second);
            for (typename vector<block>::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
            {
                r(get<0>(*it), get<1>(*it)) = get<2>(*it);
            }
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_left(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(1, rcd.second);
            for (typename vector<block>::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
            {
                r(0, get<1>(*it)) = get<2>(*it);
            }
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_right(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(rcd.first, 1);
            for (typename vector<block>::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
            {
                r(get<0>(*it), 0) = get<2>(*it);
            }
            return r;
        }
    };

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
                      , boost::shared_ptr<OPTable<Matrix, SymmGroup> > op_tags_)
        : length(length_)
        , used_dims(length)
        , ident(ident_)
        , tag_prempo(length)
        , maximum(2)
        , finalized(false)
        , leftmost_right(length)
        , op_tags(op_tags_)
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

                // retrieve the actual operator from the tag table
                op_t current_op = (*op_tags)[it->second];

                // apply scaling factor to first operator
                scale_type cur_scale = 1.;
                if (it == ops.begin()) {
                    current_op *= term.scale;
                    cur_scale = term.scale;
                }

                if (minp != maxp) { // bond term
                    /** tag_prempo filling **/
                    tag_prempo[it->first].push_back(boost::make_tuple(first_use_b, second_use_b, it->second, cur_scale));

                    used_dims[it->first].insert(use_b);
                } else // site term
                    site_terms[it->first] += current_op;
                done[it->first] = true;
            }
            
            // put fill ops in between the nontrivial ops  (only needed for 1- and 2-terms)
            if (term.operators.size() <= 2) {
                op_t fill = (*op_tags)[term.fill_operator];
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

        boost::shared_ptr<OPTable<Matrix, SymmGroup> > op_tags;
        vector<vector<tag_block> > tag_prempo;
        
        pos_t maximum, leftmost_right;

        void finalize()
        {
            for (typename std::map<pos_t, op_t>::const_iterator it = site_terms.begin();
                 it != site_terms.end(); ++it) {

                // TODO implement plus operation
                tag_type site_tag = op_tags->register_site_op(it->second);
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
            MPOTensor<Matrix, SymmGroup> r(rcd.first, rcd.second, ops, op_tags);
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_left(vector<tag_block> const & ops)
        {
            pair<std::size_t, std::size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(1, rcd.second, ops, op_tags);
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_right(vector<tag_block> const & ops)
        {
            pair<std::size_t, std::size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(rcd.first, 1, ops, op_tags);
            return r;
        }
    };
    
    
    template<class Matrix, class SymmGroup>
    class CorrMakerBase {
    public:
        virtual MPO<Matrix, SymmGroup> create_mpo()=0;
        virtual std::string description () const=0;
        virtual vector<vector<size_t> > numeric_labels()=0;
    };

    template<class Matrix, class SymmGroup>
    class CorrMaker : public CorrMakerBase<Matrix, SymmGroup>
    {
        typedef block_matrix<Matrix, SymmGroup> op_t;
        typedef boost::tuple<size_t, size_t, op_t> block;
        typedef vector<
        pair<
        block_matrix<Matrix, SymmGroup>,
        block_matrix<Matrix, SymmGroup>
        >
        > op_pairs;
        typedef boost::tuple<size_t, size_t, string> tag;
        
    public:
        CorrMaker(std::size_t L,
                  op_t const & identity_,
                  op_t const & fill_,
                  std::vector<std::pair<op_t, bool> > const & ops_,
                  int ref = -1)
        : prempo(L)
        , tags(L)
        , used(L)
        , with_sign(L+2)
        , identity(identity_)
        , fill(fill_)
        , ops(ops_)
        {
            with_sign[0][0] = false;
        	recurse(0, 0, 0, vector<size_t>(), ref);
        }
        
        MPO<Matrix, SymmGroup> create_mpo()
        {
            MPO<Matrix, SymmGroup> r(prempo.size());
            for (size_t p = 1; p < prempo.size()-1; ++p)
                r[p] = as_bulk(prempo[p]);
            r[0] = as_left(prempo[0]);
            r[prempo.size()-1] = as_right(*prempo.rbegin());
            
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
        
        vector<vector<size_t> > numeric_labels() { return labels; }
        
    private:
        vector<vector<block> > prempo;
        vector<vector<tag> > tags;
        vector<vector<size_t> > labels;
        
        vector<set<size_t> > used;
        vector<map<size_t, bool> > with_sign;
        
        op_t identity, fill;
        vector<std::pair<op_t, bool> > ops;
        
        size_t term(size_t p, size_t u1, std::pair<op_t, bool> const & op_p, bool trivial)
        {
            std::string lab;
            op_t op;
            if (trivial) {
            	op = (with_sign[p][u1]) ? fill : identity;
            	lab = (with_sign[p][u1]) ? "filling" : "ident";
            } else {
				lab = "nontriv";
            	if (!with_sign[p][u1] && op_p.second) {
					gemm(fill, op_p.first, op);
					lab += "*fill";
				} else if (with_sign[p][u1] && !op_p.second) {
					gemm(fill, op_p.first, op);
					lab += "*fill";
				} else {
					op = op_p.first;
				}
            }
            
        	size_t u2 = 0;
            while (used[p].count(u2) > 0) ++u2;
            prempo[p].push_back( boost::make_tuple(u1, u2, op) );
            used[p].insert(u2);
           	with_sign[p+1][u2] = (op_p.second) ? !with_sign[p][u1] : with_sign[p][u1];
            //            maquis::cout << "Adding a " << lab << " term at " << p << ", " << u1 << " -> " << u2 << std::endl;
            //            maquis::cout << op;
            if (trivial)
                tags[p].push_back( boost::make_tuple(u1, u2, lab) );
            else
                tags[p].push_back( boost::make_tuple(u1, u2, lab) );
            return u2;
        }
        
        void recurse(size_t p0, size_t which, size_t use, vector<size_t> label, int ref)
        {
            if (p0 + ops.size() - which < prempo.size()) {
                size_t use_next = term(p0, use, std::make_pair(identity, false), true);
                recurse(p0+1, which, use_next, label, ref);
            }
            
            {
                if (ref >= 0 && which == 0 && p0 != ref)
                    return;
                
                size_t use_next = term(p0, use, ops[which], false);
                
                vector<size_t> label_(label);
                label_.push_back(p0);
                
                if (which == ops.size()-1) {
                    size_t t1 = use_next, t2 = use_next;
                    for (size_t p2 = p0+1; p2 < prempo.size(); ++p2) {
                        t2 = term(p2, t1, std::make_pair(identity, false), true);
                        t1 = t2;
                    }
                    labels.resize(std::max(t2+1, labels.size()));
                    labels[t2] = label_;
                } else {
                    recurse(p0+1, which+1, use_next, label_, ref);
                }
            }
        }
        
        MPOTensor<Matrix, SymmGroup> as_bulk(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(rcd.first, rcd.second);
            for (typename vector<block>::const_iterator it = ops.begin(); it != ops.end(); ++it)
                r(get<0>(*it), get<1>(*it)) = get<2>(*it);
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_left(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(1, rcd.second);
            for (typename vector<block>::const_iterator it = ops.begin(); it != ops.end(); ++it)
                r(0, get<1>(*it)) = get<2>(*it);
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_right(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(rcd.first, rcd.second);
            for (typename vector<block>::const_iterator it = ops.begin(); it != ops.end(); ++it)
                r(get<0>(*it), get<1>(*it)) = get<2>(*it);
            return r;
        }
    };
    
    // same as CorrMaker, but operators in ops have to be even,
    //  and are avaluated as ops[0](i)*ops[1](i+1)*ops[2](j)*ops[3](j+1)
    template<class Matrix, class SymmGroup>
    class CorrMakerNN : public CorrMakerBase<Matrix, SymmGroup>
    {
        typedef block_matrix<Matrix, SymmGroup> op_t;
        typedef boost::tuple<size_t, size_t, op_t> block;
        typedef vector<
        pair<
        block_matrix<Matrix, SymmGroup>,
        block_matrix<Matrix, SymmGroup>
        >
        > op_pairs;
        typedef boost::tuple<size_t, size_t, string> tag;
        
    public:
        CorrMakerNN(std::size_t L,
                    op_t const & identity_,
                    op_t const & fill_,
                    std::vector<std::pair<op_t, bool> > const & ops_,
                    int ref = -1)
        : prempo(L)
        , tags(L)
        , used(L)
        , with_sign(L+2)
        , identity(identity_)
        , fill(fill_)
        , ops(ops_)
        {
            assert(ops.size() % 2 == 0);
            with_sign[0][0] = false;
            recurse(0, 0, 0, vector<size_t>(), ref);
        }
        
        MPO<Matrix, SymmGroup> create_mpo()
        {
            MPO<Matrix, SymmGroup> r(prempo.size());
            for (size_t p = 1; p < prempo.size()-1; ++p)
                r[p] = as_bulk(prempo[p]);
            r[0] = as_left(prempo[0]);
            r[prempo.size()-1] = as_right(*prempo.rbegin());
            
            return r;
        }
        
        vector<vector<size_t> > numeric_labels() { return labels; }
        
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
        
    private:
        vector<vector<block> > prempo;
        vector<vector<tag> > tags;
        vector<vector<size_t> > labels;
        
        vector<map<size_t, bool> > with_sign;
        vector<set<size_t> > used;
        
        op_t identity, fill;
        vector<std::pair<op_t, bool> > ops;
        
        size_t term(size_t p, size_t u1, std::pair<op_t, bool> const & op_p, bool trivial)
        {
            std::string lab;
            op_t op;
            if (trivial) {
            	op = (with_sign[p][u1]) ? fill : identity;
            	lab = (with_sign[p][u1]) ? "filling" : "ident";
            } else {
				lab = "nontriv";
            	if (!with_sign[p][u1] && op_p.second) {
					gemm(fill, op_p.first, op);
					lab += "*fill";
				} else if (with_sign[p][u1] && !op_p.second) {
					gemm(fill, op_p.first, op);
					lab += "*fill";
				} else {
					op = op_p.first;
				}
            }
            
        	size_t u2 = 0;
            while (used[p].count(u2) > 0) ++u2;
            prempo[p].push_back( boost::make_tuple(u1, u2, op) );
            used[p].insert(u2);
           	with_sign[p+1][u2] = (op_p.second) ? !with_sign[p][u1] : with_sign[p][u1];
            //            maquis::cout << "Adding a " << lab << " term at " << p << ", " << u1 << " -> " << u2 << std::endl;
            //            maquis::cout << op;
            if (trivial)
                tags[p].push_back( boost::make_tuple(u1, u2, lab) );
            else
                tags[p].push_back( boost::make_tuple(u1, u2, lab) );
            return u2;
        }
        
        void recurse(size_t p0, size_t which, size_t use, vector<size_t> label, int ref)
        {
            if (p0 + ops.size() - which < prempo.size()) {
                size_t use_next = term(p0, use, std::make_pair(identity, false),  true);
                recurse(p0+1, which, use_next, label, ref);
            }
            
            {
                if (ref >= 0 && which == 0 && p0 != ref)
                    return;
                
                size_t use_next = term(p0++, use, ops[which++], false);
                use_next = term(p0, use_next, ops[which], false);
                
                vector<size_t> label_(label);
                label_.push_back(p0-1);
                label_.push_back(p0);
                
                if (which == ops.size()-1) {
                    size_t t1 = use_next, t2 = use_next;
                    for (size_t p2 = p0+1; p2 < prempo.size(); ++p2) {
                        t2 = term(p2, t1, std::make_pair(identity, false), true);
                        t1 = t2;
                    }
                    labels.resize(std::max(t2+1, labels.size()));
                    labels[t2] = label_;
                } else {
                    recurse(p0+1, which+1, use_next, label_, ref);
                }
            }
        }
        
        MPOTensor<Matrix, SymmGroup> as_bulk(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(rcd.first, rcd.second);
            for (typename vector<block>::const_iterator it = ops.begin(); it != ops.end(); ++it)
                r(get<0>(*it), get<1>(*it)) = get<2>(*it);
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_left(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(1, rcd.second);
            for (typename vector<block>::const_iterator it = ops.begin(); it != ops.end(); ++it)
                r(0, get<1>(*it)) = get<2>(*it);
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_right(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(rcd.first, rcd.second);
            for (typename vector<block>::const_iterator it = ops.begin(); it != ops.end(); ++it)
                r(get<0>(*it), get<1>(*it)) = get<2>(*it);
            return r;
        }
    };
    
}


#endif
