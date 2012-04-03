// These are Hamiltonians that describe the interactions of our model
// They are sparse matrices implemented as an v = M*v operation.

#include <limits>

namespace hp2c
{

typedef std::size_t index_type;

template <typename Graph>
class sparse_matrix
{
    protected:
//        typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_type;
        typedef unsigned int vertex_type;
        typedef typename boost::graph_traits<Graph>::edge_descriptor edge_type;
        typedef typename boost::graph_traits<Graph>::edge_iterator edge_iterator;
        //
        // data members
        //
        const Graph graph_; 
    public:
        sparse_matrix(Graph const& graph)
            :graph_(graph)
        {
        }
        
        /**
          * Returns the dimension of the Hilbert space on
          * which the (square) sparse matrix acts.
          */
        std::size_t get_dimension() const
        {
            return pow(2.0,static_cast<int>(num_vertices(graph_)));
        }
        
        int symmetry_sectors_begin()
        {
            return 0;
        }

        int symmetry_sectors_end()
        {
            return 1;
        }
};


template <typename Graph>
class ising_hamiltonian : public sparse_matrix<Graph>
{
    private:
        typedef unsigned int vertex_type;
        typedef typename boost::graph_traits<Graph>::edge_iterator edge_iterator;
    public:
        ising_hamiltonian(Graph const& graph, int symmetry_sector)
            :sparse_matrix<Graph>(graph)
        {
            assert(symmetry_sector == 0);
        }

        /**
          * Does a matrix vector multiplication with the vector of polynomials.
          */
        template <typename PolynomialVector>
        PolynomialVector apply(PolynomialVector const& v)
        {

            // H =    -J * \sum_{<i,j>} \sigma^z_i \sigma^z_j     - h \sum_i \sigma^x_i
            PolynomialVector result(v.size());
            for(typename PolynomialVector::size_type  index = 0; index < v.size(); ++index)
            {
                // TODO check if v[index] == 0
                // if(v[index] == 0)
                // continue;

                // J term
                for(std::pair<edge_iterator,edge_iterator> e_r = edges(sparse_matrix<Graph>::graph_); e_r.first != e_r.second; ++e_r.first)
                {
                    vertex_type vtx_A = source(*e_r.first,sparse_matrix<Graph>::graph_);
                    vertex_type vtx_B = target(*e_r.first,sparse_matrix<Graph>::graph_);
                    unsigned int vtx_A_state = get_vertex_state_from_index(index,vtx_A);
                    unsigned int vtx_B_state = get_vertex_state_from_index(index,vtx_B);

                    if(vtx_A_state == vtx_B_state)
                    {
                        result[index] += v[index]*monomial_type(1,0); // <=> result[index] = v[index];
                    }
                    else
                    {
                        result[index] += -v[index]*monomial_type(1,0);
                    }
                }

                // h term (external field)
                for(unsigned int vtx=0; vtx < num_vertices(sparse_matrix<Graph>::graph_); ++vtx)
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
    protected:
        /**
          * Extracts the physical state of a single vertex/site from
          * the index of the vector.
          * (The state of all vertices/sites is encoded in the index.)
          */
        static inline unsigned int get_vertex_state_from_index(index_type index, vertex_type vertex)
        {
            return (index >> vertex) & index_type(1);
        }

        /**
          * Encodes the new states of vertices/sites 'vtx_A' and 'vtx_B'
          * in the index 'index'.
          * (The information about all other vertices in 'index' is kept.)
          */
        static inline index_type generate_index_from_index_and_vertex_state_pair (
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
        static inline index_type generate_index_from_index_and_vertex_state (
                index_type index,
                vertex_type vtx,
                unsigned int vtx_state
                )
        {
            index_type bitmask = (index_type(1) << vtx);
            index_type changes = index_type(vtx_state) << vtx;
            return (index & ~bitmask) | changes;
        }

};


template <typename Graph>
class heisenberg_hamiltonian : public sparse_matrix<Graph>
{
    private:
        typedef unsigned int vertex_type;
        typedef typename boost::graph_traits<Graph>::edge_iterator edge_iterator;
    public:
        heisenberg_hamiltonian(Graph const& graph, int symmetry_sector)
            : sparse_matrix<Graph>(graph), index_(1<<num_vertices(graph))
        {
            // Generate basis
            unsigned int basis_size=0;
            for(unsigned int i=0; i < index_.size(); ++i)
            {
                const int mag = magnetization(i,num_vertices(graph));
                if(mag == symmetry_sector)
                {
                    states_.push_back(i);
                    index_[i] = states_.size() - 1;
                    ++basis_size;
                }
                else
                {
                    index_[i] = std::numeric_limits<index_type>::max();
                }
            }
        }

        int symmetry_sectors_begin()
        {
            return -static_cast<int>(num_vertices(sparse_matrix<Graph>::graph_));
        }

        int symmetry_sectors_end()
        {
            return static_cast<int>(num_vertices(sparse_matrix<Graph>::graph_))+1;
        }

        int magnetization(unsigned int state, unsigned int num_vertices) const
        {
            int bitcount = 0;
            for(unsigned int i=0; i < num_vertices; ++i)
            {
                bitcount += (state >> i) & 1u;
            }
            return (2*bitcount)-static_cast<int>(num_vertices);
        }

        std::size_t get_dimension() const
        {
            return states_.size();
        }

        template <typename PolynomialVector>
        PolynomialVector apply(PolynomialVector const& v)
        {
            PolynomialVector result(v.size());
            for(typename PolynomialVector::size_type index=0; index <v.size(); ++index)
            {
                // Get the state in a bit representation
                const unsigned int state = get_state_from_index(index);

                for(std::pair<edge_iterator,edge_iterator> e_r = edges(sparse_matrix<Graph>::graph_); e_r.first != e_r.second; ++e_r.first)
                {
                    vertex_type vtx_A = source(*e_r.first,sparse_matrix<Graph>::graph_);
                    vertex_type vtx_B = target(*e_r.first,sparse_matrix<Graph>::graph_);
                    const unsigned int vtx_A_state = get_vertex_state_from_state(state,vtx_A);
                    const unsigned int vtx_B_state = get_vertex_state_from_state(state,vtx_B);
                    assert( (vtx_A_state & ~1u) == 0);
                    assert( (vtx_B_state & ~1u) == 0);

                    if(vtx_A_state != vtx_B_state)
                    {
                        // Sz Sz
                        result[index] += -v[index]*monomial_type(1,0);

                        // S+ S-
                        result[generate_index_from_state_and_vertex_state_pair(state, vtx_A, vtx_B, vtx_B_state, vtx_A_state)] += 2*v[index]*monomial_type(1,0);
                    }
                    else
                    {
                        result[index] += v[index];
                    }
                }
            }
            return result;
        }
    private:
        std::vector<index_type> index_;
        std::vector<unsigned int> states_;

        /**
          * Get state from index
          */
        inline unsigned int get_state_from_index(index_type i)
        {
            assert(i < states_.size());
            return states_[i];
        }

        /**
          * Extracts the physical state of a single vertex/site from
          * the index of the vector.
          * (The state of all vertices/sites is encoded in the index.)
          */
        static inline unsigned int get_vertex_state_from_state(unsigned int state, vertex_type vertex)
        {
            return (state >> vertex) & 1u;
        }

        /**
          * Encodes the new states of vertices/sites 'vtx_A' and 'vtx_B'
          * in the index 'index'.
          * (The information about all other vertices in 'index' is kept.)
          */
        inline index_type generate_index_from_state_and_vertex_state_pair (
                unsigned int old_state,
                vertex_type vtx_A,
                vertex_type vtx_B,
                unsigned int vtx_A_state,
                unsigned int vtx_B_state
                )
        {
            unsigned int bitmask = (1u << vtx_A) | (1u << vtx_B);
            unsigned int changes = (vtx_A_state << vtx_A) | (vtx_B_state << vtx_B);
            return index_[(old_state & ~bitmask) | changes];
        }
       
        /**
          * Same as generate_index_from_index_and_vertex_state_pair()
          * but for a single vertex.
          */ 
        inline index_type generate_index_from_index_and_vertex_state (
                unsigned int old_state,
                vertex_type vtx,
                unsigned int vtx_state
                )
        {
            unsigned int bitmask = (1u << vtx);
            unsigned int changes = vtx_state << vtx;
            return index_[(old_state & ~bitmask) | changes];
        }
};
} // namespace hp2c
