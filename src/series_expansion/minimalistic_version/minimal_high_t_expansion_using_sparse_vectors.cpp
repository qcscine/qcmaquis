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
 *                       calculated the program will generate a chain of
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
#include "minimal_polynomial.hpp"

// A large integer of 128-256 bits (fixed size)
// this will define type >  large_int  <
#ifdef USE_GMP_INTEGERS
#include "integer_classes/use_gmp_integers.hpp"
#endif //USE_GMP_INTEGERS

#ifdef USE_INT64_INTEGERS
namespace hp2c
{
    typedef int64_t large_int;
}
#endif //USE_INT64_INTEGERS




using namespace hp2c;

//
// Typedefs
//


// A polynomial with large_int coefficients
typedef polynomial<large_int> polynomial_type;

// A sparse vector of polynomials (filling about 10%-15%)
typedef std::size_t index_type;
typedef std::map<index_type,polynomial_type> sparse_vector_type;





/**
  * Returns the index of the element of the sparse vector to which iterator 'it' points.
  */
inline index_type get_index_from_iterator(sparse_vector_type::const_iterator const& it)
{
    return it->first;
}



/**
  * The inner product of two sparse vectors. (We assume same dimension.)
  */
polynomial_type inner_product(sparse_vector_type const& a, sparse_vector_type const& b)
{
    if(a.size() > b.size())
    {
        return inner_product(b,a);
    }

    polynomial_type result;
    for(sparse_vector_type::const_iterator it = a.begin(); it != a.end(); ++it)
    {
        index_type index = get_index_from_iterator(it);
        sparse_vector_type::const_iterator it_b = b.find(index);
        if(it_b != b.end())
            result += it->second * it_b->second;
    }
    return result;
}


class sparse_matrix
{
    private:
        // Information for the sparse matrix (given at runtime)
        typedef unsigned int vertex_type;
        typedef std::pair<vertex_type, vertex_type> edge_type;
        
        const std::vector<edge_type> edge_list;
        const unsigned int num_vertices;
        

        /**
          * Extracts the physical state of a single vertex/site from
          * the index of the sparse vector.
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

    public:
        sparse_matrix(std::vector<edge_type> const& edge_list, unsigned int num_vertices)
            :edge_list(edge_list),num_vertices(num_vertices)
        {
        }

        /**
          * Returns the dimension of the Hilbert space on
          * which the (square) sparse matrix acts.
          */
        std::size_t get_dimension() const
        {
            return std::pow(2,static_cast<unsigned int>(num_vertices));
        }

        /**
          * Does a matrix vector multiplication with a sparse vector.
          */
        sparse_vector_type apply(sparse_vector_type const& v)
        {

            // H =    -J * \sum_{<i,j>} \sigma^z_i \sigma^z_j     - h \sum_i \sigma^x_i
            sparse_vector_type result;
            for(sparse_vector_type::const_iterator it = v.begin(); it != v.end(); ++it)
            {
                index_type index = get_index_from_iterator(it);

                // J term
                for(
                    std::vector<edge_type>::const_iterator edge_it
                        = edge_list.begin();
                    edge_it != edge_list.end();
                    ++edge_it
                    )
                {
                    vertex_type vtx_A = edge_it->first;
                    vertex_type vtx_B = edge_it->second;
                    unsigned int vtx_A_state = get_vertex_state_from_index(index,vtx_A);
                    unsigned int vtx_B_state = get_vertex_state_from_index(index,vtx_B);

                    if(vtx_A_state == vtx_B_state)
                    {
                        result[index] += it->second*monomial<large_int>(1,0); // <=> result[index] = v[index];
                    }
                    else
                    {
                        result[index] += -1*monomial<large_int>(1,0)*it->second;
                    }
                }

                // h term (external field)
                for(unsigned int vtx=0; vtx < num_vertices; ++vtx)
                {
                    // it->second == v[index]
                    unsigned int vtx_state = get_vertex_state_from_index(index,vtx);
                    if( vtx_state == 0 )
                    {
                        result[
                            generate_index_from_index_and_vertex_state(index,vtx,1)
                            ] += it->second*monomial<large_int>(0,1);
                    }
                    else
                    {
                        assert(vtx_state == 1);
                        result[
                            generate_index_from_index_and_vertex_state(index,vtx,0)
                            ] += it->second*monomial<large_int>(0,1);
                    }

                } 
            }

            return result;
        }
};

class high_t_expansion
{
    private:
        // The truncation order of the series expansion
        const unsigned int max_order;
        
        // Information for the sparse matrix
        typedef unsigned int vertex_type;
        typedef std::pair<vertex_type, vertex_type> edge_type;

        const unsigned int num_vertices;
        std::vector<edge_type> edge_list;
    public:

        high_t_expansion(unsigned int max_order,unsigned int num_vertices)
            : max_order(max_order), num_vertices(num_vertices)
        {
            // generate a simple graph as example
            for(vertex_type vtx = 1; vtx < num_vertices; ++vtx)
                edge_list.push_back(edge_type(vtx-1,vtx));
        }

        polynomial_type exec()
        {
            sparse_matrix sp_matrix(edge_list,num_vertices);

            polynomial_type result;
            result += sp_matrix.get_dimension(); // 0th order contribution
            for(std::size_t i = 0; i < sp_matrix.get_dimension(); ++i)
            {
                // Create an initial sparse vector with a single entry
                sparse_vector_type state;
                state[i] = (polynomial_type()+= 1);

                for(unsigned int n=1; n < max_order/2; ++n)
                {
                    sparse_vector_type previous_state(state);

                    state = sp_matrix.apply(state);

                    result += inner_product(state,state);
                    result += inner_product(state,previous_state);
                }
            }
            return result;
        }

};

large_int faculty(int i)
{
    if (i <= 1)
        return 1;
    return i*faculty(i-1);
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
    for(std::size_t j = 0; j < r.max_order; ++j)
    {
        for(std::size_t h = 0; h < r.max_order; ++h)
        {
            if(r(j,h) > 0)
                std::cout <<" +"<< r(j,h) << "/"<< faculty(j+h)<<"*J^"<<j<<"*h^"<<h;
            else if(r(j,h) < 0)
                std::cout <<" "<< r(j,h) << "/"<< faculty(j+h)<<"*J^"<<j<<"*h^"<<h;
        }
    }

    std::cout<<std::endl;
    return 0;
}
