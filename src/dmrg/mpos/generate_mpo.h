#ifndef GENERATE_MPO_H
#define GENERATE_MPO_H

#include "adjancency.h"
#include "mp_tensors/mpo.h"

#include "block_matrix/block_matrix.h"
#include "block_matrix/block_matrix_algorithms.h"

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
        
        virtual void push_extra_terms(MPOMaker<Matrix, SymmGroup>&, Adjacency&) { }
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
        MPOMaker(Adjacency const & adj_,
                 Hamiltonian<Matrix, SymmGroup> & H_)
        : adj(adj_)
        , H(H_)
        , used_dims(adj.size())
        , prempo(adj.size())
        , maximum(0)
        {   
            for (size_t p = 0; p < adj.size(); ++p)
            {
                if (p+1 < adj.size())
                    prempo[p].push_back( make_tuple(0, 0, H.get_identity()) );
                if (p > 0)
                    prempo[p].push_back( make_tuple(1, 1, H.get_identity()) );
            }
        }
        
        void add_bond_ops()
        {
            for (size_t p = 0; p < adj.size(); ++p) {
                vector<size_t> neighs = adj[p];
                for (vector<size_t>::iterator neigh = neighs.begin(); neigh != neighs.end(); ++neigh)
                    for (int i = 0; i < H.num_2site_ops(); ++i)
                    {
                        pair<op_t, op_t> op = H.get_2site_op(i);
                        
                        vector<pair<size_t, op_t> > ops;
                        ops.push_back( make_pair( p, op.first ) );
                        ops.push_back( make_pair( *neigh, op.second ) );
                        add_term(ops);
                    }
            }
            
            cout << "Maximum: " << maximum << endl;
        }
        
        void add_term(vector<pair<size_t, op_t> > & ops)
        {
            std::sort(ops.begin(), ops.end(), compare<Matrix, SymmGroup>);
            
            vector<size_t> positions;
            for (typename vector<pair<size_t, op_t> >::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
                positions.push_back( it->first );
//            std::copy(positions.begin(), positions.end(),
//                      std::ostream_iterator<size_t>(cout, " ")); cout << endl;
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

            vector<bool> done(adj.size(), false);
            for (typename vector<pair<size_t, op_t> >::const_iterator it = ops.begin();
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
                    prempo[p].push_back( make_tuple(use_b, use_b, H.get_free()) );
                    used_dims[p].insert(use_b);
                    done[p] = true;
                }
        }
            
        MPO<Matrix, SymmGroup> create_mpo()
        {
            MPO<Matrix, SymmGroup> r(adj.size());
            for (size_t p = 1; p < adj.size() - 1; ++p)
                r[p] = as_bulk(prempo[p], maximum+1);
            r[0] = as_left(prempo[0], maximum+1);
            r[adj.size()-1] = as_right(prempo[adj.size()-1], maximum+1);
            
            return r;
        }
        
    private:
        Adjacency const & adj;
        Hamiltonian<Matrix, SymmGroup> & H;

        vector<set<size_t> > used_dims;
        vector<vector<block> > prempo;
        
        size_t maximum;
        
        MPOTensor<Matrix, SymmGroup> as_bulk(vector<block> const & ops,
                                             size_t M)
        {
            MPOTensor<Matrix, SymmGroup> r(M, M);
            for (typename vector<block>::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
            {
                r(get<0>(*it), get<1>(*it)) = get<2>(*it);
            }
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_left(vector<block> const & ops,
                                             size_t M)
        {
            MPOTensor<Matrix, SymmGroup> r(1, M);
            for (typename vector<block>::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
            {
                r(0, get<1>(*it)) = get<2>(*it);
            }
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_right(vector<block> const & ops,
                                              size_t M)
        {
            MPOTensor<Matrix, SymmGroup> r(M, 1);
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
