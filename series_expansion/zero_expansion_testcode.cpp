#include <boost/graph/adjacency_list.hpp>
#include "numeric/simple_sparse_vector/simple_sparse_vector.hpp"
#include "numeric/simple_polynomial.hpp"
#include <boost/rational.hpp>
#include "expansion_modules/zero_expansion.hpp"
#include "expansion_modules/high_t_expansion.hpp"
#include "operators.hpp"
#include <vector>
#include <iostream>

// For GMP
#include <gmpxx.h>
#include "util/gmp_limits.h"

#ifndef MAXORDER
#define MAXORDER 8
#endif

namespace series_expansion
{

template <typename Operator, typename Vector>
class degenerate_state_generator
{
	public:
        typedef typename Vector::size_type index_type;

        template <typename Graph>
 		std::vector<index_type> operator() (const Graph& g, unsigned int sector)
		{
			std::vector<index_type> result;
			index_type b(0);
			if(sector == 0)
			{
				result.push_back(b);
			}
			else if(sector == 1)
			{
				for(unsigned int i=0; i<num_vertices(g); ++i)
				{
					result.push_back(index_from_site_state_pair<Operator::states_per_site>(b,i,Operator::trp_minus));
				}
			}
			return result;
		}
};

}


using namespace series_expansion;

template <typename RationalNumType, typename Graph>
void do_zero_expansion(Graph const& g)
{
    zero_expansion<
       heisenberg_h0,
       heisenberg_v,
       degenerate_state_generator<heisenberg_h0,simple_sparse_vector<RationalNumType> >,
       Graph,
       simple_sparse_vector<RationalNumType>,
       simple_polynomial<RationalNumType>
           >
                zero_exp(heisenberg_h0(),heisenberg_v(),MAXORDER);
    
    blas::dense_matrix<simple_polynomial<RationalNumType> > ground_state_result = zero_exp(g,0);
    blas::dense_matrix<simple_polynomial<RationalNumType> > excited_state_result = zero_exp(g,1);
    std::cout<<ground_state_result<<std::endl<<std::endl;
    std::cout<<excited_state_result<<std::endl<<std::endl;
}

int main()
{
    //
    // What do we want to use as large integers?
    //
    typedef int large_int_type;
    //typedef mpz_class large_int_type;   // GMP

    //
    // What do we want to use as rational numbers?
    //
    //typedef boost::rational<large_int_type> rational_num_type;
    typedef mpq_class rational_num_type;  // GMP



	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph_type;

    graph_type g;

    // generate some graph with about 15 edges
    add_edge(0,1,g);
    add_edge(0,2,g);
    add_edge(0,3,g);
    add_edge(1,4,g);
    add_edge(2,5,g);
    add_edge(4,6,g);
    add_edge(5,6,g);
    add_edge(3,7,g);
    add_edge(6,9,g);
    add_edge(7,8,g);
    
    do_zero_expansion<rational_num_type>(g);

    return 0;
}

