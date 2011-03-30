#include <valarray>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <cassert>
//
// Typedefs
//
#define NOT_IN_USE false

// A large integer of 128-256 bits (fixed size)
typedef int large_int;

// A polynomial with large_int coefficients
typedef std::valarray<large_int> polynomial_type;

// A sparse vector of polynomials (filling about 10%-15%)
typedef std::size_t index_type;
typedef std::map<index_type,polynomial_type> sparse_vector_type;

inline index_type get_index_from_iterator(sparse_vector_type::const_iterator const& it)
{
    return it->first;
}

polynomial_type inner_product(sparse_vector_type const& a, sparse_vector_type const& b)
{
    if(a.size() > b.size())
    {
        return inner_product(b,a);
    }

    polynomial_type result(1,0);
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
        enum{ VERTICES = 5 };
        std::vector<edge_type> edge_list;
        

        inline unsigned int get_vertex_state_from_index(index_type index, vertex_type vertex)
        {
            return (index >> vertex) & index_type(1);
        }

        inline index_type generate_index_from_index_and_vertex_state_pair (
                index_type index,
                vertex_type vtx_A,
                vertex_type vtx_B,
                unsigned int vtx_A_state,
                unsigned int vtx_B_state
                )
        {
            assert( NOT_IN_USE );
            index_type bitmask = (index_type(1) << vtx_A) | (index_type(1) << vtx_B);
            index_type changes = index_type(vtx_A_state) << vtx_A | index_type(vtx_B_state) << vtx_B;
            return (index & ~bitmask) | changes;
        }
        
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
        sparse_matrix()
        {
            // generate a simple graph as example
            for(vertex_type vtx = 1; vtx <= VERTICES; ++vtx)
                edge_list.push_back(edge_type(vtx-1,vtx));
        }

        std::size_t get_dimension()
        {
            return std::pow(2,static_cast<unsigned int>(VERTICES));
        }

        sparse_vector_type apply(sparse_vector_type const& v)
        {

            // H =    -J * \sum_{<i,j>} \sigma^z_i \sigma^z_j     - h \sum_i \sigma^x_i
            sparse_vector_type result;
            for(sparse_vector_type::const_iterator it = v.begin(); it != v.end(); ++it)
            {
                index_type index = get_index_from_iterator(it);

                // J term
                for(
                    std::vector<edge_type>::iterator edge_it
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
                        result[index] += it->second; // <=> result[index] = v[index];
                    }
                    else
                    {
//                        index_type second_index =
//                            generate_index_from_index_and_vertex_states(
//                                    index,
//                                    vtx_A,
//                                    vtx_B,
//                                    vtx_B_state,
//                                    vtx_A_state
//                                    );
//
//                        // it->second == v[index]
//                        result[index] = -1*it->second;
//                        result[second_index] = 2*it->second;
                        result[index] += -1*it->second;
                    }
                }

                // h term (external field)
                for(unsigned int vtx=0; vtx <= VERTICES; ++vtx)
                {
                    // it->second == v[index]
                    unsigned int vtx_state = get_vertex_state_from_index(index,vtx);
                    if( vtx_state == 0 )
                    {
                        result[
                            generate_index_from_index_and_vertex_state(index,vtx,1)
                            ] += it->second;
                    }
                    else
                    {
                        assert(vtx_state == 1);
                        result[
                            generate_index_from_index_and_vertex_state(index,vtx,0)
                            ] += it->second;
                    }

                } 
            }

            // TODO finish implementation of matrix for transverse field ising
            return result;
        }
};

class high_t_expansion
{
    private:
        // The dimension of the sparse vector of polynomials
        unsigned int max_order;
    public:

        high_t_expansion(unsigned int max_order)
            : max_order(max_order)
        {
            // Create a simple bond list for the sparse matrix
            // We will have about a million of these lists.
//            sparse_vector_dim = 2;
//            for(unsigned int i=0; i < MAX_ORDER; ++i)
//            {
//                bond_list.push_back(i,i+1);
//                sparse_vector_dim *= 2;
//            }
        }

        polynomial_type exec()
        {
            sparse_matrix sp_matrix;
            polynomial_type result(max_order);
            for(std::size_t i = 0; i < sp_matrix.get_dimension(); ++i)
            {
                // Create an initial sparse vector with a single entry
                sparse_vector_type state;
                state[i] = polynomial_type(1,1);

                for(unsigned int n=1; n < max_order/2; ++n)
                {
                    sparse_vector_type previous_state(state);

                    state = sp_matrix.apply(state);

//                    result[2*n]     += inner_product(state,state);
//                    result[2*n-1]   += inner_product(state,previous_state);
                    result += inner_product(state,state);
                    result += inner_product(state,previous_state);
                }
            }
            return result;
        }

};

int main()
{
    high_t_expansion ht(5);
    polynomial_type r = ht.exec();
    for(std::size_t i = 0; i != r.size(); ++i)
        std::cout << r[i] << "+";
    std::cout<<std::endl;
    return 0;
}
