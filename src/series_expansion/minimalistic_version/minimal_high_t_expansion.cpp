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
 *          ./a.out <expansion_order> <num_vertices>
 *
 *      expansion_order: the maximal order of the series expansion.
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
#include <cstdlib>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <cassert>

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
    typedef hp2c::polynomial<large_int> polynomial_type;
    typedef std::vector<polynomial_type> polynomial_vector_type;
}
#endif //USE_INT64_INTEGERS

#ifdef USE_GMP_INTEGERS
#include "integer_classes/use_gmp_integers.hpp"
#include "minimal_polynomial.hpp"
namespace hp2c
{
    typedef mpz_class large_int;
    typedef hp2c::monomial<large_int> monomial_type;
    typedef hp2c::polynomial<large_int> polynomial_type;
    typedef std::vector<polynomial_type> polynomial_vector_type;
}
#endif //USE_GMP_INTEGERS

#ifdef USE_VLI_INTEGERS_CPU

#ifdef VLI_USE_GPU
#include "vli/utils/gpu_manager.h"
#include "vli/utils/gpu_manager.hpp"
#endif //VLI_USE_GPU

#include "vli/polynomial/vector_polynomial_cpu.hpp"
#include "vli/polynomial/polynomial_cpu.hpp"
#include "vli/polynomial/monomial.hpp"
#include "vli/vli_cpu.hpp"
#include "vli/vli_traits.hpp"
namespace hp2c
{
    typedef vli::vli_cpu<unsigned long int,3> large_int;
    typedef vli::monomial<large_int> monomial_type;
    typedef vli::polynomial_cpu<large_int, POLYNOMIAL_MAX_ORDER > polynomial_type;
    typedef vli::vector_polynomial_cpu< polynomial_type > polynomial_vector_type;
}
#endif //USE_VLI_INTEGERS_CPU

// ------------------------------------------------------------------
// Code
// ------------------------------------------------------------------

using namespace hp2c;

typedef std::size_t index_type;

class sparse_matrix
{
    private:
        typedef unsigned int vertex_type;
        typedef std::pair<vertex_type, vertex_type> edge_type;
    public:
        sparse_matrix(std::vector<edge_type> const& edge_list, unsigned int num_vertices)
            :edge_list_(edge_list),num_vertices_(num_vertices)
        {
        }

        /**
          * Returns the dimension of the Hilbert space on
          * which the (square) sparse matrix acts.
          */
        std::size_t get_dimension() const
        {
            return pow(2,static_cast<unsigned int>(num_vertices_));
        }

        /**
          * Does a matrix vector multiplication with the vector of polynomials.
          */
        polynomial_vector_type apply(polynomial_vector_type const& v)
        {

            // H =    -J * \sum_{<i,j>} \sigma^z_i \sigma^z_j     - h \sum_i \sigma^x_i
            polynomial_vector_type result(v.size());
            for(polynomial_vector_type::size_type  index = 0; index < v.size(); ++index)
            {
                // TODO check if v[index] == 0
                // if(v[index] == 0)
                // continue;

                // J term
                for(
                    std::vector<edge_type>::const_iterator edge_it
                        = edge_list_.begin();
                    edge_it != edge_list_.end();
                    ++edge_it
                    )
                {
                    vertex_type vtx_A = edge_it->first;
                    vertex_type vtx_B = edge_it->second;
                    unsigned int vtx_A_state = get_vertex_state_from_index(index,vtx_A);
                    unsigned int vtx_B_state = get_vertex_state_from_index(index,vtx_B);

                    if(vtx_A_state == vtx_B_state)
                    {
                        result[index] += v[index]*monomial_type(1,0); // <=> result[index] = v[index];
                    }
                    else
                    {
                        result[index] += -1*(monomial_type(1,0)*v[index]); //negate result will be very faster
                    }
                }

                // h term (external field)
                for(unsigned int vtx=0; vtx < num_vertices_; ++vtx)
                {
                    // it->second == v[index]
                    unsigned int vtx_state = get_vertex_state_from_index(index,vtx);
                    if( vtx_state == 0 )
                    {
                        index_type r = generate_index_from_index_and_vertex_state(index,vtx,1);
                        result[r] += v[index]*monomial_type(0,1);
                    }
                    else
                    {
                        //assert(vtx_state == 1);
                        index_type r = generate_index_from_index_and_vertex_state(index,vtx,0);
                        result[r] += v[index]*monomial_type(0,1);
                    }

                }
            }

            return result;
        }
    private:
        /**
          * Extracts the physical state of a single vertex/site from
          * the index of the vector.
          * (The state of all vertices/sites is encoded in the index.)
          */
        inline unsigned int get_vertex_state_from_index(index_type index, vertex_type vertex)
        {
            return (index >> vertex) & index_type(1);
        }

        /**
          * Encodes the new states of vertices/sites 'vtx_A' and 'vtx_B'
          * in the index 'index'.
          * (The information about all other vertices in 'index' is kept.)
          */
        inline index_type generate_index_from_index_and_vertex_state_pair (
                index_type index,
                vertex_type vtx_A,
                vertex_type vtx_B,
                unsigned int vtx_A_state,
                unsigned int vtx_B_state
                )
        {
            index_type bitmask = (index_type(1) << vtx_A) | (index_type(1) << vtx_B);
            index_type changes = index_type(vtx_A_state) << vtx_A | index_type(vtx_B_state) << vtx_B;
            return (index & ~bitmask) | changes;
        }
       
        /**
          * Same as generate_index_from_index_and_vertex_state_pair()
          * but for a single vertex.
          */ 
        inline index_type generate_index_from_index_and_vertex_state (
                index_type index,
                vertex_type vtx,
                unsigned int vtx_state
                )
        {
            index_type bitmask = (index_type(1) << vtx);
            index_type changes = index_type(vtx_state) << vtx;
            return (index & ~bitmask) | changes;
        }

        // 
        // Information for the sparse matrix (given at runtime)
        //
        
        const std::vector<edge_type> edge_list_;
        const unsigned int num_vertices_;
        

};

class high_t_expansion
{
    private:
        // Information for the sparse matrix
        typedef unsigned int vertex_type;
        typedef std::pair<vertex_type, vertex_type> edge_type;
    
    public:
        high_t_expansion(unsigned int max_order,unsigned int num_vertices)
            : max_order_(max_order), num_vertices_(num_vertices)
        {
            // generate a simple graph as example
            for(vertex_type vtx = 1; vtx < num_vertices_; ++vtx)
                edge_list_.push_back(edge_type(vtx-1,vtx));
        }

        polynomial_type exec()
        {
            sparse_matrix sp_matrix(edge_list_,num_vertices_);

            polynomial_type result;
            result += sp_matrix.get_dimension(); // 0th order contribution
            // Compute the trace of (sparse_matrix)^max_order
            for(std::size_t i = 0; i < sp_matrix.get_dimension(); ++i)
            {
                // Create an initial vector with a single entry
                polynomial_vector_type state(sp_matrix.get_dimension());
                state[i] = (polynomial_type()+= 1);
                for(unsigned int n=1; n < max_order_/2; ++n)
                {
                    polynomial_vector_type previous_state(state);

                    state = sp_matrix.apply(state);

                    result += inner_product(state,state);
                    result += inner_product(state,previous_state);
                }
            }
            return result;
        }

    private:
        // The truncation order of the series expansion
        const unsigned int max_order_;
        const unsigned int num_vertices_;
        std::vector<edge_type> edge_list_;
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
    //
    // Runtime parameters
    //
    unsigned int expansion_order = 20;
    unsigned int num_vertices = 8;
    if(argc == 3)
    {
        // Default values may be overwritten by command line arguments.
        expansion_order = atoi(argv[1]);
        num_vertices = atoi(argv[2]);
    }

    
    
    // Run expansion
    high_t_expansion ht(expansion_order,num_vertices);
    polynomial_type r = ht.exec();
    
    
    
    // Print out the result
    for(unsigned int j = 0; j < r.max_order; ++j)
    {
        for(unsigned int h = 0; h < r.max_order; ++h)
        {
            if(r(j,h) > 0)
                std::cout <<" +"<< r(j,h) << "/"<< factorial(j+h)<<"*J^"<<j<<"*h^"<<h;
            else if(r(j,h) < 0)
                std::cout <<" "<< r(j,h) << "/"<< factorial(j+h)<<"*J^"<<j<<"*h^"<<h << std::endl;
        }
    }

    std::cout<<std::endl;
    return 0;
}
