/*****************************************************************************
 *
 * A minimal implementation of a high-temperature series expansion of a
 *
 * Ising chain with transverse field 'h'.
 *
 * (C) 2010 by Andreas Hehn (hehn@phys.ethz.ch)
 *
 * The program calculates the contribution of ONE cluster to the partition
 * function.
 *
 * Parameters:
 *  Compile time (macro):
 *    POLYNOMIAL_MAX_ORDER : the order at which polynomials are truncated
 *
 *  Runtime:
 *
 *          ./a.out <num_vertices>
 *
 *      num_vertices:    the number of vertices/sites of the cluster to be
 *                       calculated. the program will generate a chain of
 *                       'num_vertices' sites.
 *          
 *
 * In a real application we will have a few million clusters(=graphs) with
 * about 20 vertices and 20 edges each. The structure of the graphs will
 * be more complicated than the chain presented here, but this won't not change
 * the computation too much.
 * We would like to calculate the contribution of each cluster up to at least
 * 20-th order for this model.
 * 
*
 *****************************************************************************/
#include <boost/graph/adjacency_list.hpp>
#include <cstdlib>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <cassert>
#include <string>
#include "utils/timings.h"
// ------------------------------------------------------------------
// Typedefs and includes for the different integer classes
// ------------------------------------------------------------------
//
// A large integer of 128-256 bits (fixed size)
// this will define type >  large_int  <

#ifdef USE_INT64_INTEGERS
namespace hp2c
{
    typedef int64_t large_int;
    typedef hp2c::monomial<int> monomial_type;
}
#endif //USE_INT64_INTEGERS

#ifdef USE_GMP_INTEGERS
#include "integer_classes/use_gmp_integers.hpp"
#include "minimal_polynomial.hpp"
namespace hp2c
{
    typedef mpz_class large_int;
    typedef mpz_class double_large_int;
    typedef hp2c::monomial<large_int> monomial_type;
}
#endif //USE_GMP_INTEGERS

#ifdef USE_LUKAS_VLI_INTEGERS
#include <alps/graph/vli.hpp>
#include "minimal_polynomial.hpp"
namespace hp2c
{
    typedef alps::graph::vli<256> large_int;
    typedef hp2c::monomial<large_int> monomial_type;
}
#endif //USE_LUKAS_VLI_INTERGERS


#ifdef USE_VLI_INTEGERS_CPU

#ifdef VLI_USE_GPU
#include "vli/utils/gpu_manager.h"
#include "vli/utils/gpu_manager.hpp"
#endif //VLI_USE_GPU

#include "vli/polynomial/vector_polynomial_cpu.hpp"
#include "vli/polynomial/polynomial_cpu.h"
#include "vli/polynomial/monomial.hpp"
#include "vli/vli_cpu.h"
#include "vli/vli_traits.hpp"
namespace hp2c
{
    typedef vli::vli_cpu<unsigned long int,3> large_int;
    typedef vli::vli_cpu<unsigned long int,6> double_large_int;
    typedef vli::monomial<large_int> monomial_type;
}
#endif //USE_VLI_INTEGERS_CPU

// ------------------------------------------------------------------
// Code
// ------------------------------------------------------------------

#include "hamiltonians.hpp"

using namespace hp2c;


template <typename Polynomial>
class polynomial_vector
{
};

#ifdef USE_VLI_INTEGERS_CPU
template <typename T, unsigned int Order>
class polynomial_vector<vli::polynomial_cpu<T,Order> > : public vli::vector_polynomial_cpu<vli::polynomial_cpu<T,Order> >
{
    public:
        polynomial_vector(std::size_t size)
            : vli::vector_polynomial_cpu<vli::polynomial_cpu<T,Order> >(size)
        {
        }
};
#else //USE_VLI_INTEGERS_CPU
template <typename T, unsigned int Order>
class polynomial_vector<hp2c::polynomial<T,Order> > : public std::vector<hp2c::polynomial<T,Order> >
{
    public:
        polynomial_vector(std::size_t size)
            : std::vector<hp2c::polynomial<T,Order> >(size)
        {
        }
};

#endif //USE_VLI_INTEGERS_CPU

template <typename Graph, typename SparseMatrix, unsigned int Order>
class high_t_expansion
{
    public:
        high_t_expansion(Graph const& g)
            : graph_(g)
        {
        }

        template <template <class,unsigned int> class Polynomial>
        Polynomial<double_large_int,Order> exec()
        {
            assert(Order % 2 == 0); //TODO compile time assert
            SparseMatrix hamiltonian(graph_,0);

            Polynomial<double_large_int,Order> result;
            // Compute the trace of (sparse_matrix)^max_order
            for(int symmetry_sector = hamiltonian.symmetry_sectors_begin(); symmetry_sector != hamiltonian.symmetry_sectors_end(); ++symmetry_sector)
            { 
//                std::cout<<"m_z="<<symmetry_sector<<"\t";
                Polynomial<double_large_int,Order> contrib;
                SparseMatrix sp_matrix(graph_,symmetry_sector);
                contrib += sp_matrix.get_dimension(); // 0th order contribution
                for(std::size_t i = 0; i < sp_matrix.get_dimension(); ++i)
                {
                    // Create an initial vector with a single entry
                    polynomial_vector<Polynomial<large_int,Order/2> > state(sp_matrix.get_dimension());
                    state[i] = (Polynomial<large_int,Order/2>() += 1);
                    for(unsigned int n=1; n < Order/2; ++n)
                    {
                        polynomial_vector<Polynomial<large_int,Order/2> > previous_state(state);

                        state = sp_matrix.apply(state);

                        contrib += inner_product(state,state);
                        contrib += inner_product(state,previous_state);
                    }
                }
//                std::cout<<contrib<<std::endl;
                result += contrib;
            }
            return result;
        }

    private:
        // The truncation order of the series expansion
        const Graph graph_;
};

large_int factorial(unsigned int i)
{
    if (i <= 1)
        return large_int(1);

    large_int a(i);
    return a*factorial(i-1);
}

int main(int argc, char** argv)
{

    // TODO
    static const unsigned int expansion_order = 20;

#ifdef VLI_USE_GPU
    gpu::gpu_manager* GPU;
    GPU->instance();
#endif


//    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::property<alps::vertex_type_t,alps::type_type>, boost::property<alps::edge_type_t,alps::type_type>  > graph_type;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph_type;

#ifdef USE_VLI_INTEGERS_CPU
    typedef vli::polynomial_cpu<large_int, expansion_order> polynomial_type;
#else
    typedef hp2c::polynomial<large_int, expansion_order> polynomial_type;
#endif
   
    //
    // Runtime parameters
    //
    unsigned int num_vertices = 8;
    std::string timer_name("default");
    if(argc == 3)
    {
        // Default values may be overwritten by command line arguments.
        num_vertices = atoi(argv[1]);
        timer_name = std::string(argv[2]);
    }
    else
    {
        std::cout<<"Using default values..."<<std::endl;
    }
    // Build graph

    graph_type g;
    add_vertex(g);
    for(unsigned int i=1; i < num_vertices; ++i)
    {
        add_vertex(g);
        add_edge(i-1,i,g);
    }
    
    Timer A(timer_name);    
    A.begin();
    
    // Run expansion
    high_t_expansion<graph_type, heisenberg_hamiltonian<graph_type>, expansion_order> ht(g);
//    high_t_expansion<graph_type, ising_hamiltonian<graph_type>, expansion_order> ht(g);
#ifdef USE_VLI_INTEGERS_CPU
    vli::polynomial_cpu<double_large_int, expansion_order> r = ht.exec<vli::polynomial_cpu>();
#else
    hp2c::polynomial<large_int, expansion_order> r = ht.exec<hp2c::polynomial>();
#endif

    A.end();
   // A.save(); 

    // Print out the result
    std::cout<<r<<std::endl;

    return 0;
}
