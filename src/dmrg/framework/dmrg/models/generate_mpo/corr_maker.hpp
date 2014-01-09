/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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

#ifndef GENERATE_MPO_CORR_MAKER_H
#define GENERATE_MPO_CORR_MAKER_H

#include "dmrg/models/generate_mpo/utils.hpp"

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"

#include "dmrg/mp_tensors/mpo.h"

#include "dmrg/models/lattice.h"
#include "dmrg/models/op_handler.h"

#include <string>
#include <sstream>

namespace generate_mpo
{
    template<class Matrix, class SymmGroup>
    class CorrMakerBase {
    public:
        virtual ~CorrMakerBase() {}
        virtual MPO<Matrix, SymmGroup> create_mpo()=0;
        virtual std::string description () const=0;
        virtual vector<vector<size_t> > const& numeric_labels()=0;
    };

    template<class Matrix, class SymmGroup>
    class CorrMaker : public CorrMakerBase<Matrix, SymmGroup>
    {
        typedef tag_detail::tag_type tag_type;
        typedef block_matrix<Matrix, SymmGroup> op_t;
        typedef boost::tuple<size_t, size_t, tag_type, typename Matrix::value_type> block;
        typedef boost::tuple<size_t, size_t, string> tag;
        
    public:
        CorrMaker(Lattice const& lat_,
                  const std::vector<op_t> & ident_,
                  const std::vector<op_t> & fill_,
                  std::vector<std::pair<std::vector<op_t>, bool> > const & ops,
                  int ref = -1)
        : lat(lat_)
        , prempo(lat.size())
        , tags(lat.size())
        , used(lat.size())
        , with_sign(lat.size()+2)
        , identities(ident_.size())
        , fillings(fill_.size())
        , op_tags(ops.size())
        {
            /// register operators
            for (int type=0; type<ident_.size(); ++type)
                identities[type] = tag_handler.register_op(ident_[type], tag_detail::bosonic);
            for (int type=0; type<fill_.size(); ++type)
                fillings[type] = tag_handler.register_op(fill_[type], tag_detail::bosonic);
            
            for (size_t n=0; n<ops.size(); ++n) {
                op_tags[n].first.resize(ops[n].first.size());
                op_tags[n].second = ops[n].second;
                std::vector<tag_type> & tops = op_tags[n].first;
                for (int type=0; type<tops.size(); ++type)
                    tops[type] = tag_handler.register_op(ops[n].first[type], ops[n].second ? tag_detail::fermionic : tag_detail::bosonic);
            }

            with_sign[0][0] = false;
        	recurse(0, 0, 0, std::vector<size_t>(), ref);
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
        Lattice const& lat;
        TagHandler<Matrix, SymmGroup> tag_handler;

        vector<vector<block> > prempo;
        vector<vector<tag> > tags;
        vector<vector<size_t> > labels;
        
        vector<set<size_t> > used;
        vector<map<size_t, bool> > with_sign;
        std::vector<tag_type> identities, fillings;
        // TODO: use just vector<tag_type>, as there is the is_fermionic() function in TagHandler
        vector<std::pair<std::vector<tag_type>, bool> > op_tags;
        
        size_t term(size_t p, size_t u1, std::pair<std::vector<tag_type>, bool> const & op_p, bool trivial)
        {
            std::string lab;
            tag_type op;
            typename Matrix::value_type scale;
            if (trivial) {
            	op = (with_sign[p][u1]) ? fillings[lat.get_prop<int>("type", p)] : identities[lat.get_prop<int>("type", p)];
                scale = 1.;
            	lab = (with_sign[p][u1]) ? "filling" : "ident";
            } else {
				lab = "nontriv";
            	if (!with_sign[p][u1] && op_p.second) {
					//gemm(fill, op_p.first, op);
					boost::tie(op, scale) = tag_handler.get_product_tag(fillings[lat.get_prop<int>("type", p)], op_p.first[lat.get_prop<int>("type", p)]);
					lab += "*fill";
				} else if (with_sign[p][u1] && !op_p.second) {
					//gemm(fill, op_p.first, op);
					boost::tie(op, scale) = tag_handler.get_product_tag(fillings[lat.get_prop<int>("type", p)], op_p.first[lat.get_prop<int>("type", p)]);
					lab += "*fill";
				} else {
					op = op_p.first[lat.get_prop<int>("type", p)];
                    scale = 1.;
				}
            }
            
        	size_t u2 = 0;
            while (used[p].count(u2) > 0) ++u2;
            prempo[p].push_back( boost::make_tuple(u1, u2, op, scale) );
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
            if (p0 + op_tags.size() - which < prempo.size()) {
                size_t use_next = term(p0, use, std::make_pair(identities, false), true);
                recurse(p0+1, which, use_next, label, ref);
            }
            
            {
                if (ref >= 0 && which == 0 && p0 != ref)
                    return;
                
                if (tag_handler.get_op(op_tags[which].first[lat.get_prop<int>("type", p0)]).n_blocks() == 0)
                    return;
                
                size_t use_next = term(p0, use, op_tags[which], false);
                
                vector<size_t> label_(label);
                label_.push_back(p0);
                
                if (which == op_tags.size()-1) {
                    size_t t1 = use_next, t2 = use_next;
                    for (size_t p2 = p0+1; p2 < prempo.size(); ++p2) {
                        t2 = term(p2, t1, std::make_pair(identities, false), true);
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
            MPOTensor<Matrix, SymmGroup> r(rcd.first, rcd.second, ops, tag_handler.get_operator_table());
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_left(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(1, rcd.second, ops, tag_handler.get_operator_table());
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_right(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(rcd.first, rcd.second, ops, tag_handler.get_operator_table());
            return r;
        }
    };
    
    // same as CorrMaker, but operators in ops have to be even,
    //  and are avaluated as ops[0](i)*ops[1](i+1)*ops[2](j)*ops[3](j+1)
    template<class Matrix, class SymmGroup>
    class CorrMakerNN : public CorrMakerBase<Matrix, SymmGroup>
    {
        typedef tag_detail::tag_type tag_type;
        typedef block_matrix<Matrix, SymmGroup> op_t;
        typedef boost::tuple<size_t, size_t, tag_type, typename Matrix::value_type> block;
        typedef boost::tuple<size_t, size_t, string> tag;

    public:
        CorrMakerNN(Lattice const& lat_,
                    const std::vector<op_t> & ident_,
                    const std::vector<op_t> & fill_,
                    std::vector<std::pair<std::vector<op_t>, bool> > const & ops,
                    int ref = -1)
        : lat(lat_)
        , prempo(lat.size())
        , tags(lat.size())
        , used(lat.size())
        , with_sign(lat.size()+2)
        , identities(ident_.size())
        , fillings(fill_.size())
        , op_tags(ops.size())
        {
            assert(ops.size() % 2 == 0);

            /// register operators
            for (int type=0; type<ident_.size(); ++type)
                identities[type] = tag_handler.register_op(ident_[type], tag_detail::bosonic);
            for (int type=0; type<fill_.size(); ++type)
                fillings[type] = tag_handler.register_op(fill_[type], tag_detail::bosonic);

            for (size_t n=0; n<ops.size(); ++n) {
                op_tags[n].first.resize(ops[n].first.size());
                op_tags[n].second = ops[n].second;
                std::vector<tag_type> & tops = op_tags[n].first;
                for (int type=0; type<tops.size(); ++type)
                    tops[type] = tag_handler.register_op(ops[n].first[type], ops[n].second ? tag_detail::fermionic : tag_detail::bosonic);
            }

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
        Lattice const& lat;
        TagHandler<Matrix, SymmGroup> tag_handler;

        vector<vector<block> > prempo;
        vector<vector<tag> > tags;
        vector<vector<size_t> > labels;
        
        vector<set<size_t> > used;
        vector<map<size_t, bool> > with_sign;
        
        std::vector<tag_type> identities, fillings;
        // TODO: use just vector<tag_type>, as there is the is_fermionic() function in TagHandler
        vector<std::pair<std::vector<tag_type>, bool> > op_tags;

        size_t term(size_t p, size_t u1, std::pair<std::vector<tag_type>, bool> const & op_p, bool trivial)
        {
            std::string lab;
            tag_type op;
            typename Matrix::value_type scale;
            if (trivial) {
            	op = (with_sign[p][u1]) ? fillings[lat.get_prop<int>("type", p)] : identities[lat.get_prop<int>("type", p)];
            	lab = (with_sign[p][u1]) ? "filling" : "ident";
            } else {
				lab = "nontriv";
            	if (!with_sign[p][u1] && op_p.second) {
					//gemm(fill, op_p.first, op);
					boost::tie(op, scale) = tag_handler.get_product_tag(fillings[lat.get_prop<int>("type", p)], op_p.first[lat.get_prop<int>("type", p)]);
					lab += "*fill";
				} else if (with_sign[p][u1] && !op_p.second) {
					//gemm(fill, op_p.first, op);
					boost::tie(op, scale) = tag_handler.get_product_tag(fillings[lat.get_prop<int>("type", p)], op_p.first[lat.get_prop<int>("type", p)]);
					lab += "*fill";
				} else {
					op = op_p.first[lat.get_prop<int>("type", p)];
				}
            }
            
        	size_t u2 = 0;
            while (used[p].count(u2) > 0) ++u2;
            prempo[p].push_back( boost::make_tuple(u1, u2, op, 1.0) );
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
            if (p0 + op_tags.size() - which < prempo.size()) {
                size_t use_next = term(p0, use, std::make_pair(identities, false),  true);
                recurse(p0+1, which, use_next, label, ref);
            }
            
            {
                if (ref >= 0 && which == 0 && p0 != ref)
                    return;
                
                if (tag_handler.get_op(op_tags[which].first[lat.get_prop<int>("type", p0)]).n_blocks() == 0)
                    return;
                size_t use_next = term(p0++, use, op_tags[which++], false);

                if (tag_handler.get_op(op_tags[which].first[lat.get_prop<int>("type", p0)]).n_blocks() == 0)
                    return;
                use_next = term(p0, use_next, op_tags[which], false);
                
                vector<size_t> label_(label);
                label_.push_back(p0-1);
                label_.push_back(p0);
                
                if (which == op_tags.size()-1) {
                    size_t t1 = use_next, t2 = use_next;
                    for (size_t p2 = p0+1; p2 < prempo.size(); ++p2) {
                        t2 = term(p2, t1, std::make_pair(identities, false), true);
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
            MPOTensor<Matrix, SymmGroup> r(rcd.first, rcd.second, ops, tag_handler.get_operator_table());
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_left(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(1, rcd.second, ops, tag_handler.get_operator_table());
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_right(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(rcd.first, rcd.second, ops, tag_handler.get_operator_table());
            return r;
        }
    };
    
}


#endif
