/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef GENERATE_MPO_H
#define GENERATE_MPO_H

#include "block_matrix/block_matrix.h"
#include "block_matrix/block_matrix_algorithms.h"
#include "block_matrix/symmetry.h"

#include "mp_tensors/mpo.h"

#include "lattice.h"

#include <string>
#include <sstream>

namespace generate_mpo
{
    
	template<class Matrix, class SymmGroup>
	struct Operator_Term
	{
		typedef block_matrix<Matrix, SymmGroup> op_t;
        
		std::vector<std::pair<typename Lattice::pos_t, op_t> > operators;
		op_t fill_operator;
        
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
    
    template<class Matrix, class SymmGroup>
    void compress_on_bond(std::vector<boost::tuple<size_t, size_t, block_matrix<Matrix, SymmGroup> > > & pm1,
                          std::vector<boost::tuple<size_t, size_t, block_matrix<Matrix, SymmGroup> > > & pm2)
    {
        using namespace std;
        using namespace boost::tuples;
        typedef block_matrix<Matrix, SymmGroup> op_t;
        typedef boost::tuple<size_t, size_t, op_t> block;
        
        set<size_t> bond_used_dims;
        for (typename vector<block>::iterator it = pm1.begin(); it != pm1.end(); ++it)
            bond_used_dims.insert(get<1>(*it));
        for (typename vector<block>::iterator it = pm2.begin(); it != pm2.end(); ++it)
            bond_used_dims.insert(get<0>(*it));
        
        //        cout << "Compression: " << *max_element(bond_used_dims.begin(),
        //                                                bond_used_dims.end()) << " -> " << bond_used_dims.size() << endl;
        
        map<size_t, size_t> compression_map;
        size_t c = 0;
        for (set<size_t>::iterator it = bond_used_dims.begin();
             it != bond_used_dims.end(); ++it)
            compression_map[*it] = c++;
        
        for (typename vector<block>::iterator it = pm1.begin(); it != pm1.end(); ++it)
            get<1>(*it) = compression_map[get<1>(*it)];
        for (typename vector<block>::iterator it = pm2.begin(); it != pm2.end(); ++it)
            get<0>(*it) = compression_map[get<0>(*it)];
    }
    
    template<class Matrix, class SymmGroup>
    std::pair<size_t, size_t>
    rcdim(std::vector<boost::tuple<size_t, size_t, block_matrix<Matrix, SymmGroup> > > const & pm)
    {
        using namespace std;
        using namespace boost::tuples;
        typedef block_matrix<Matrix, SymmGroup> op_t;
        typedef boost::tuple<size_t, size_t, op_t> block;
        
        list<size_t> l, r;
        
        for (typename vector<block>::const_iterator it = pm.begin(); it != pm.end(); ++it) {
            l.push_back( get<0>(*it) );
            r.push_back( get<1>(*it) );
        }
        
        return make_pair(*max_element(l.begin(), l.end())+1,
                         *max_element(r.begin(), r.end())+1);
    }
    
    template<class Matrix, class SymmGroup>
    bool compare(pair<size_t, block_matrix<Matrix, SymmGroup> > const & p1,
                 pair<size_t, block_matrix<Matrix, SymmGroup> > const & p2)
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
        , maximum(0)
        , leftmost_right(length)
        {   
            for (size_t p = 0; p < length; ++p)
            {
                if (p+1 < length)
                    prempo[p].push_back( make_tuple(0, 0, ident) );
            }
        }
        
        void add_term(Operator_Term<Matrix, SymmGroup> const & term)
        {
            // TODO: removed const&, because of sorting (non-const operation)
            std::vector<std::pair<typename Lattice::pos_t, op_t> > ops = term.operators;
            
            std::sort(ops.begin(), ops.end(), compare<Matrix, SymmGroup>);
            
            vector<size_t> positions;
            for (typename vector<pair<typename Lattice::pos_t, op_t> >::const_iterator
                 it = ops.begin();
                 it != ops.end(); ++it)
                positions.push_back( it->first );
            size_t minp = *min_element(positions.begin(), positions.end());
            size_t maxp = *max_element(positions.begin(), positions.end());
            
            size_t use_b;
            for (use_b = 2; true; ++use_b) {
                bool valid = true;
                for (size_t p = minp; p < maxp; ++p)
                    valid &= (used_dims[p].count(use_b) == 0);
                if (valid)
                    break;
            }
            maximum = max(use_b, maximum);
            
            vector<bool> done(length, false);
            for (typename vector<pair<typename Lattice::pos_t, op_t> >::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
            {
                size_t first_use_b = (it->first == minp ? 0 : use_b);
                size_t second_use_b = (it->first == maxp ? 1 : use_b);
                assert( it->first < prempo.size() );
                prempo[it->first].push_back( make_tuple(first_use_b, second_use_b, it->second) );
                if (minp != maxp) { // bond term
                    prempo[it->first].push_back( make_tuple(first_use_b, second_use_b, it->second) );
                    used_dims[it->first].insert(use_b);
                } else // site term
                    site_terms[it->first] += it->second;
                done[it->first] = true;
            }
            
            for (size_t p = minp; p <= maxp; ++p)
                if (!done[p]) {
                    prempo[p].push_back( make_tuple(use_b, use_b, term.fill_operator) );
                    used_dims[p].insert(use_b);
                    done[p] = true;
                }
            
            leftmost_right = std::min(leftmost_right, maxp);
        }
        
        MPO<Matrix, SymmGroup> create_mpo()
        {
            for (typename std::map<std::size_t, op_t>::const_iterator it = site_terms.begin();
                 it != site_terms.end(); ++it)
                prempo[it->first].push_back( make_tuple(0, 1, it->second) );
            
            for (size_t p = leftmost_right + 1; p < length; ++p)
                prempo[p].push_back( make_tuple(1, 1, ident) );
//            
//            for (typename vector<vector<block> >::iterator it = prempo.begin();
//                 it + 1 != prempo.end();
//                 ++it)
//                compress_on_bond(*it, *(it+1));
            
            MPO<Matrix, SymmGroup> r(length);
            for (size_t p = 1; p < length - 1; ++p)
                r[p] = as_bulk(prempo[p]);
            r[0] = as_left(prempo[0]);
            r[length-1] = as_right(prempo[length-1]);
            
            return r;
        }
        
    private:
        std::size_t length;
        block_matrix<Matrix, SymmGroup> ident;
        vector<set<size_t> > used_dims;
        vector<vector<block> > prempo;
        std::map<std::size_t, op_t> site_terms;
        
        size_t maximum, leftmost_right;
        
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
    class CorrMaker
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
            prempo[p].push_back( make_tuple(u1, u2, op) );
            used[p].insert(u2);
           	with_sign[p+1][u2] = (op_p.second) ? !with_sign[p][u1] : with_sign[p][u1];
            //            cout << "Adding a " << lab << " term at " << p << ", " << u1 << " -> " << u2 << endl;
            //            cout << op;
            if (trivial)
                tags[p].push_back( make_tuple(u1, u2, lab) );
            else
                tags[p].push_back( make_tuple(u1, u2, lab) );
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
    class CorrMakerNN
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
            prempo[p].push_back( make_tuple(u1, u2, op) );
            used[p].insert(u2);
           	with_sign[p+1][u2] = (op_p.second) ? !with_sign[p][u1] : with_sign[p][u1];
            //            cout << "Adding a " << lab << " term at " << p << ", " << u1 << " -> " << u2 << endl;
            //            cout << op;
            if (trivial)
                tags[p].push_back( make_tuple(u1, u2, lab) );
            else
                tags[p].push_back( make_tuple(u1, u2, lab) );
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
