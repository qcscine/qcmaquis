#include <boost/graph/adjacency_list.hpp>
#include "numeric/simple_sparse_vector/simple_sparse_vector.hpp"
#include "numeric/simple_polynomial.hpp"
#include "expansion_modules/zero_expansion.hpp"
#include "expansion_modules/high_t_expansion.hpp"
#include "operators.hpp"
#include <vector>
#include <iostream>

// For GMP use
#include <gmpxx.h>
#include "util/gmp_limits.h"


#ifndef MAXORDER
#define MAXORDER 10
#endif

using namespace series_expansion;


template <typename LargeIntType, typename Graph>
void do_high_t_expansion(Graph const& g)
{
    high_t_expansion<
        heisenberg_hamiltonian,
        Graph,
        simple_sparse_vector<LargeIntType>,
        simple_polynomial<LargeIntType>
            >
            high_t_exp(heisenberg_hamiltonian(),MAXORDER);
    std::pair<simple_polynomial<LargeIntType>,std::vector<LargeIntType> > result = high_t_exp(g);

    // Result output
    for(unsigned int i =0; i < result.first.size(); ++i)
    {
        if(result.first[i] >= 0)
            std::cout<<"+";
        std::cout<<result.first[i]<<"/"<<result.second[i]<<"*l^"<<i;
    }
    std::cout<<std::endl;
}

int main()
{
    //
    // What do we want to use as large integers?
    //
    //typedef int large_int_type;
    typedef mpz_class large_int_type;   // GMP
	
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph_type;

    graph_type g;

    // generate some graph with about 15 edges
    add_edge(0,1,g);
    add_edge(0,2,g);
    add_edge(0,3,g);
    add_edge(1,4,g);
    add_edge(2,5,g);
//    add_edge(4,6,g);
//    add_edge(5,6,g);
//    add_edge(3,7,g);
//    add_edge(6,9,g);
//    add_edge(7,8,g);
    
    do_high_t_expansion<large_int_type>(g);

    return 0;
}

