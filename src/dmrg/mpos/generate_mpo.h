#ifndef GENERATE_MPO_H
#define GENERATE_MPO_H

#include "adjancency.h"
#include "mp_tensors/mpo.h"

#include <boost/tuple/tuple.hpp>

#include <set>

namespace mpos {
    using namespace std;
    using namespace boost::tuples;
    
    size_t next_free(std::vector<std::size_t> const & out_taken,
                     std::vector<std::size_t> const & in_taken)
    {
        for (std::size_t k = 0; true; ++k)
        {
            if (count(out_taken.begin(), out_taken.end(), k) == 0 &&
                count(in_taken.begin(), in_taken.end(), k) == 0)
                return k;
        }
    }
    
    size_t next_free(set<size_t> const & s)
    {
        for (std::size_t k = 2; true; ++k)
            if (s.count(k) == 0)
                return k;
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
        static op_pairs hb_ops(double Jxy, double Jz)
        {
            block_matrix<Matrix, U1> ident, splus, sminus, sz, zero;
            
            ident.insert_block( make_tuple(Matrix(1, 1, 1), -1, -1) );
            ident.insert_block( make_tuple(Matrix(1, 1, 1), 1, 1) );
            
            splus.insert_block( make_tuple(Matrix(1, 1, 1), -1, 1) );
            
            sminus.insert_block( make_tuple(Matrix(1, 1, 1), 1, -1) );
            
            sz.insert_block( make_tuple(Matrix(1, 1, 0.5), 1, 1) );
            sz.insert_block( make_tuple(Matrix(1, 1, -0.5), -1, -1) );
            
            op_pairs ret;
            
            ret.push_back(make_pair(ident, ident));
            ret.push_back(make_pair(Jxy/2*splus, sminus));
            ret.push_back(make_pair(Jxy/2*sminus, splus));
            ret.push_back(make_pair(Jz*sz, sz));
            
            return ret;
        }
        
        static MPOTensor<Matrix, SymmGroup> as_bulk(vector<block> const & ops,
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
        
        static MPOTensor<Matrix, SymmGroup> as_left(vector<block> const & ops,
                                                    size_t M)
        {
            MPOTensor<Matrix, SymmGroup> r(1, M);
            for (typename vector<block>::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
            {
                if (get<1>(*it) == 1)
                    continue;
                r(0, get<1>(*it)) = get<2>(*it);
            }
            return r;
        }
        
        static MPOTensor<Matrix, SymmGroup> as_right(vector<block> const & ops,
                                                     size_t M)
        {
            MPOTensor<Matrix, SymmGroup> r(M, 1);
            for (typename vector<block>::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
            {
                if (get<0>(*it) == 0)
                    continue;
                r(get<0>(*it), 0) = get<2>(*it);
            }
            return r;
        }
        
        static MPO<Matrix, SymmGroup>
        create_mpo(Adjacency const & adj,
                   op_pairs const & ops)
        {
            vector<set<size_t> > used_dims(adj.size()-1);
            
            vector<vector<block> > prempo(adj.size());
            size_t maximum = 0;
            
            for (size_t p = 0; p < adj.size(); ++p)
            {
                prempo[p].push_back( make_tuple(0, 0, ops[0].first) );
                prempo[p].push_back( make_tuple(1, 1, ops[0].first) );
                
                vector<size_t> neighs = adj[p];
                for (vector<size_t>::iterator neigh = neighs.begin(); neigh != neighs.end(); ++neigh)
                    for (typename op_pairs::const_iterator op = ops.begin()+1;
                         op != ops.end(); ++op)
                    {
                        size_t nfree = next_free(used_dims[p]);
                        maximum = std::max(maximum, nfree);
                        
                        prempo[p].push_back( make_tuple(0, nfree, op->first) );
                        used_dims[p].insert(nfree);
                        
                        for (size_t pp = p+1; pp < *neigh; ++pp) {
                            prempo[pp].push_back( make_tuple(nfree, nfree, ops[0].first) );
                            used_dims[pp].insert(nfree);
                        }
                        
                        prempo[*neigh].push_back( make_tuple(nfree, 1, op->second) );
                    }
            }
            
            cout << "Maximum: " << maximum << endl;
            
            MPO<Matrix, SymmGroup> r(adj.size());
            for (size_t p = 1; p < adj.size() - 1; ++p)
                r[p] = as_bulk(prempo[p], maximum+1);
            r[0] = as_left(prempo[0], maximum+1);
            r[adj.size()-1] = as_right(prempo[adj.size()-1], maximum+1);
            
            return r;
        }
    };
}

#endif
