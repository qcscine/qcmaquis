/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef GENERATE_MPO_CORR_MAKER_H
#define GENERATE_MPO_CORR_MAKER_H

#include "dmrg/models/generate_mpo/utils.hpp"

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"

#include "dmrg/mp_tensors/mpo.h"

#include "dmrg/models/lattice.h"

#include <string>
#include <sstream>

namespace generate_mpo
{
    template<class Matrix, class SymmGroup>
    class CorrMakerBase {
    public:
        virtual MPO<Matrix, SymmGroup> create_mpo()=0;
        virtual std::string description () const=0;
        virtual vector<vector<size_t> > const& numeric_labels()=0;
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
        
        vector<vector<size_t> > const& numeric_labels() { return labels; }
        
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
        
        vector<vector<size_t> > const& numeric_labels() { return labels; }
        
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
