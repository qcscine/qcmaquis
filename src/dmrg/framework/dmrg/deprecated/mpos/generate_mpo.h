/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef GENERATE_MPO_H
#define GENERATE_MPO_H

#include "adjacency.h"
#include "dmrg/mp_tensors/mpo.h"

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

#include <boost/tuple/tuple.hpp>
#include <set>

namespace mpos {
    template<class Matrix, class SymmGroup> class MPOMaker;
    
    template<class Matrix, class SymmGroup>
    class Hamiltonian
    {
    public:
        typedef block_matrix<Matrix, SymmGroup> op_t;
        typedef std::pair<op_t, op_t> op_pair;
        
        virtual ~Hamiltonian() { };
        
        virtual op_t get_identity() = 0;
        virtual op_t get_free() = 0;
        
        virtual Index<SymmGroup> get_phys() = 0;
        
        virtual int num_2site_ops() = 0;
        virtual op_pair get_2site_op(int) = 0;
        
        virtual int num_1site_ops() = 0;
        virtual op_t get_1site_op(int) = 0;
        
        virtual void push_extra_terms(MPOMaker<Matrix, SymmGroup>&, adj::Adjacency&) { }
    };
    
    using namespace std;
    using namespace boost::tuples;
    
    size_t next_free(vector<size_t> const & out_taken,
                     vector<size_t> const & in_taken)
    {
        for (size_t k = 0; true; ++k)
        {
            if (count(out_taken.begin(), out_taken.end(), k) == 0 &&
                count(in_taken.begin(), in_taken.end(), k) == 0)
                return k;
        }
    }
    
    size_t next_free(set<size_t> const & s)
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
                  std::vector<op_t> const & ops_,
                  int ref = -1)
        : prempo(L)
        , tags(L)
        , used(L)
        , identity(identity_)
        , fill(fill_)
        , ops(ops_)
        {
            recurse(0, 0, 0, vector<size_t>(), ref);
        }
        
        MPO<Matrix, SymmGroup> create_mpo()
        {
//            for (typename vector<vector<block> >::iterator it = prempo.begin();
//                 it + 1 != prempo.end();
//                 ++it)
//                compress_on_bond(*it, *(it+1));
            
//            for (vector<vector<size_t> >::iterator it = labels.begin();
//                 it != labels.end(); ++it)
//            {
//                cout << "Correlator ";
//                std::copy(it->begin(), it->end(), std::ostream_iterator<size_t>(cout, " "));
//                cout << " at " << it-labels.begin() << endl;
//            }
            
//            for (size_t p = 0; p < prempo.size(); ++p)
//            {
//                cout << "Site: " << p << endl;
//                for (typename vector<tag>::const_iterator it = tags[p].begin(); it != tags[p].end(); ++it)
//                    cout << "    " << get<0>(*it) << " " << get<1>(*it) << " " << get<2>(*it) << endl;
//            }
            
            MPO<Matrix, SymmGroup> r(prempo.size());
            for (size_t p = 1; p < prempo.size()-1; ++p)
                r[p] = as_bulk(prempo[p]);
            r[0] = as_left(prempo[0]);
            r[prempo.size()-1] = as_right(*prempo.rbegin());
            
            return r;
        }
        
        vector<string> label_strings()
        {
            vector<string> ret(labels.size());
            vector<string>::iterator ot = ret.begin();
            for (vector<vector<size_t> >::iterator it = labels.begin();
                 it != labels.end(); ++it)
            {
                ostringstream oss;
                oss << "( ";
                for (vector<size_t>::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
                    oss << *it2;
                    if (it2 + 1 != it->end())
                        oss << " ) -- ( ";
                }
                oss << " )";
                *(ot++) = oss.str();
            }
            assert( ot == ret.end() );
            return ret;
        }
        
        vector<vector<size_t> > numeric_labels() { return labels; }
        
    private:
        vector<vector<block> > prempo;
        vector<vector<tag> > tags;
        vector<vector<size_t> > labels;
        
        vector<set<size_t> > used;
        
        op_t identity, fill;
        vector<op_t> ops;
        
        size_t term(size_t p, size_t u1, op_t const & op, bool trivial)
        {
            size_t u2 = 0;
            while (used[p].count(u2) > 0) ++u2;
            prempo[p].push_back( make_tuple(u1, u2, op) );
            used[p].insert(u2);
//            cout << "Adding a " << (trivial ? "trivial " : "") << "term at " << p << ", " << u1 << " -> " << u2 << endl;
            if (trivial)
                tags[p].push_back( make_tuple(u1, u2, "trivial") );
            else
                tags[p].push_back( make_tuple(u1, u2, "nontriv") );
            return u2;
        }
        
        void recurse(size_t p0, size_t which, size_t use, vector<size_t> label, int ref)
        {
            if (p0 + ops.size() - which < prempo.size()) {
                size_t use_next = term(p0, use,
                                       which == 0 ? identity : fill,
                                       true);
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
                        t2 = term(p2, t1, identity, true);
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
        MPOMaker(adj::Adjacency const & adj_,
                 Hamiltonian<Matrix, SymmGroup> & H_)
        : adj(adj_)
        , H(H_)
        , used_dims(adj.size())
        , prempo(adj.size())
        , maximum(0)
        , leftmost_right(adj.size())
        {   
            for (size_t p = 0; p < adj.size(); ++p)
            {
                if (p+1 < adj.size())
                    prempo[p].push_back( make_tuple(0, 0, H.get_identity()) );
            }
        }
        
        void add_bond_ops()
        {
            for (int p = 0; p < adj.size(); ++p) {
                vector<int> neighs = adj.forward(p);
                for (vector<int>::iterator neigh = neighs.begin(); neigh != neighs.end(); ++neigh)
                    for (int i = 0; i < H.num_2site_ops(); ++i)
                    {
                        pair<op_t, op_t> op = H.get_2site_op(i);
                        
                        vector<pair<int, op_t> > ops;
                        ops.push_back( make_pair( p, op.first ) );
                        ops.push_back( make_pair( *neigh, op.second ) );
                        add_term(ops);
                    }
            }
        }
        
        void add_term(vector<pair<int, op_t> > & ops)
        {
            std::sort(ops.begin(), ops.end(), compare<Matrix, SymmGroup>);
            
            vector<int> positions;
            for (typename vector<pair<int, op_t> >::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
                positions.push_back( it->first );
//            std::copy(positions.begin(), positions.end(),
//                      std::ostream_iterator<size_t>(cout, " ")); cout << " " << endl;
            int minp = *min_element(positions.begin(), positions.end());
            int maxp = *max_element(positions.begin(), positions.end());
            
            int use_b;
            for (use_b = 2; true; ++use_b) {
                bool valid = true;
                for (int p = minp; p < maxp; ++p)
                    valid &= (used_dims[p].count(use_b) == 0);
                if (valid)
                    break;
            }
            maximum = max(use_b, maximum);

            vector<bool> done(adj.size(), false);
            for (typename vector<pair<int, op_t> >::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
            {
                size_t first_use_b = (it->first == minp ? 0 : use_b);
                size_t second_use_b = (it->first == maxp ? 1 : use_b);
                assert( it->first < prempo.size() );
                prempo[it->first].push_back( make_tuple(first_use_b, second_use_b, it->second) );
                used_dims[it->first].insert(use_b);
                done[it->first] = true;
            }
            
            for (int p = minp; p <= maxp; ++p)
                if (!done[p]) {
                    prempo[p].push_back( make_tuple(use_b, use_b, H.get_free()) );
                    used_dims[p].insert(use_b);
                    done[p] = true;
                }
            
            leftmost_right = std::min(leftmost_right, maxp);
        }
            
        MPO<Matrix, SymmGroup> create_mpo()
        {
            for (int p = leftmost_right + 1; p < adj.size(); ++p)
                prempo[p].push_back( make_tuple(1, 1, H.get_identity()) );
            
            for (typename vector<vector<block> >::iterator it = prempo.begin();
                 it + 1 != prempo.end();
                 ++it)
                compress_on_bond(*it, *(it+1));
            
            MPO<Matrix, SymmGroup> r(adj.size());
            for (int p = 1; p < adj.size() - 1; ++p)
                r[p] = as_bulk(prempo[p]);
            r[0] = as_left(prempo[0]);
            r[adj.size()-1] = as_right(prempo[adj.size()-1]);
            
            return r;
        }
        
    private:
        adj::Adjacency const & adj;
        Hamiltonian<Matrix, SymmGroup> & H;

        vector<set<size_t> > used_dims;
        vector<vector<block> > prempo;
        
        int maximum, leftmost_right;
        
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
}

#endif
