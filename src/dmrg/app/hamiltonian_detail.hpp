/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#include "mp_tensors/mpo.h"

namespace hamiltonian_detail
{
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
        
        void add_term(app::Hamiltonian_Term<Matrix, SymmGroup> const & term)
        {
            // TODO: removed const & because of sorting (non-const operation)
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
                used_dims[it->first].insert(use_b);
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
            for (size_t p = leftmost_right + 1; p < length; ++p)
                prempo[p].push_back( make_tuple(1, 1, ident) );
            
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
}
