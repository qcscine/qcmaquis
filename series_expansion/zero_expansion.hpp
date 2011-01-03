#ifndef HP2C__ZERO_EXPANSION_HPP
#define HP2C__ZERO_EXPANSION_HPP

#include <dense_matrix/dense_matrix.h>
#include <vector>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/construct.hpp>

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

#include "numeric/simple_polynomial.hpp"
#include "numeric/inner_product.hpp"
#include "numeric/outer_product.hpp"
#include "util/calc_hilbert_space_dimension.hpp"


namespace series_expansion
{

using blas::dense_matrix;
using blas::outer_product;
#ifdef EXPANSION_DEBUG
using std::cerr;
using std::endl;
#endif //EXPANSION_DEBUG

//
//	INTERFACE
//


template <typename OperatorHo, typename OperatorV, typename InitStateGenFuncObj, typename Graph, typename Vector, typename Polynomial>
class zero_expansion
{
    private:
        OperatorHo Ho;
        OperatorV V;
        unsigned int max_order;
    public:
        typedef typename Vector::value_type value_type;

        zero_expansion(OperatorHo const& Ho, OperatorV const& V, unsigned int max_order)
            : Ho(Ho), V(V), max_order(max_order)
            {
                BOOST_STATIC_ASSERT(( boost::is_same<typename Vector::value_type, typename Polynomial::value_type>::value ));
            }
        ~zero_expansion()
            {}


        dense_matrix<Polynomial> operator () (Graph const& g, unsigned int sector) const
        {
            //
            // Vector is a vector (sparse, dense or whatever is useful)
            // to represent the quantum state |x> on the graph g
            //
            // it will contain polynoms of rational numbers (using large integers) or
            // single rational numbers as elements.
            //


            //
            // Preparation
            //
            
            #ifdef TRACE_STATE_SIZE
                // This flag enables a counter that counts how many entries the sparse vector has
                // and prints out the percentage of the full hilbert space dimension.
                std::vector<std::vector<unsigned int> > trace_state_size(max_order+1);
            #endif //TRACE_STATE_SIZE

            // Calculate the dimension the Vector must have to represent the state of the cluster.
            // (doesn't matter for the simple sparse_vector)
            typename Vector::size_type h_space_dim = calc_hilbert_space_dimension(typename Vector::size_type(),Ho,g);

            // Get a list of indices of basis states being members of the degenerate manifold
            // for this cluster and excitation number sector.
            // eg. a single excitation of a certain type
            InitStateGenFuncObj init_state_gen;
            std::vector<typename Vector::size_type> index_list(init_state_gen(g,sector));
            
            // Generate all degenerate States of the cluster from the index list
            blas::vector<Vector> init(index_list.size(),Vector(h_space_dim));
            for(unsigned int i=0; i < index_list.size(); ++i)
                init[i](index_list[i]) = value_type(1);
            
            // Calculate the eigen energy of the first of the given initial states
            value_type Eo( inner_product(init[0], Ho.apply(g,init[0]) ) );
            
            // Check if all given init states have the same eigen-energy (= are the states really degenerate?)
            for(typename std::vector<Vector>::iterator it(init.begin()+1); it != init.end(); ++it)
            {
                if(Eo != inner_product(*it, Ho.apply(g,*it)) )
                    std::runtime_error("Not all given inital states have the same energy!");
            }

            // size of the vector with the different degenerate states | init >
            // will also be the dimension of the square matrix 'Heff'
            const unsigned int num_of_degen_states = init.size();

            // Use the vector of states init now as the zero-th element of a list psi of such vectors
            // where each element psi[n] is the vector of states for the order n.
            // (Hence init is order n=0.)
            std::vector<blas::vector<Vector> > psi
                                             ( max_order+1, 
                                               blas::vector<Vector>( num_of_degen_states ,Vector(h_space_dim) )
                                             );
            psi[0] = init;
            
            // The resulting energy term Heff[order](l,m) here = Heff[order](m,l) in literature
            // Now we create a similar list for the effective Hamiltonian which we want to compute.
            // Each entry is a matrix.
            // Heff[0] == zero-th order result
            std::vector<dense_matrix<value_type> > Heff
                                                ( max_order+1,
                                                  dense_matrix<value_type>( num_of_degen_states , num_of_degen_states )
                                                );

            // The unperturbed H[0](l,m) is a diagonal matrix with the eigen-energy of the degenerate states
            for(unsigned int i=0; i<num_of_degen_states ; ++i)
                Heff[0](i,i) = Eo;

            // TODO: using tricks you only need n/2
            for(unsigned int n=1; n<=max_order; ++n)
            {
                // First recursive equation
                // move this loop further down, use vector notation
                // |psi[n][l] > = V | psi[n-1][l] > needed in all recursive equations
                for(unsigned int i=0; i < num_of_degen_states; ++i)
                    psi[n][i] = V.apply(g,psi[n-1][i]);
//                std::transform (
//                        psi[n-1].begin(),
//                        psi[n-1].end(),
//                        psi[n].begin(),
//                        boost::lambda::bind(&(V::apply), V, g, boost::lambda::_1)
//                        );
                
                // The Contributions of this order (n) are the inner products of
                // _each_ state in psi[n] with _each_ state of init.
                // So in principle its a (generalized) outer product of inner products of the states.
                Heff[n] = outer_product(psi[n],init);
                
                for(unsigned int l=0; l < num_of_degen_states; ++l)
                {
                    #ifdef EXPANSION_DEBUG
                        cout<<"--V*psi["<<n-1<<"]["<<l<<"]:--"<<endl;
                        psi[n][l].print();
                        cout<<"-------------"<<endl;
                        for(unsigned int i=0; i<L; ++i);
                        {
                            cout<<"H["<<n<<"]("<<l<<","<<i<<"):"<<H[n](l,i)<<endl;
                        }
                    #endif //EXPANSION_DEBUG

                    // TODO: replace following loops by
                    // psi[n] -=  matrix_vector_prod(Heff[m],psi[n-m]);
                    //
                    // vector of states  -= Matrix * vector of states
                    // doesn't work (for typeof(states) == dense vectors),
                    // since it is a vector of a vector and there is no
                    // generic multiplication that results in a vector
                    // of vectors, with a size promotion of the contained vectors (=states).
                    for(unsigned int m=1; m<=n; ++m)
                        for(unsigned int i=0; i<num_of_degen_states; ++i)
                            psi[n][l] -= Heff[m](l,i) * psi[n-m][i];

                    // Normalize each entry of the state with 1/(E_0 - E_k)
                    // where E_k is the eigen energy of the basis state.
                    // (= for each index of a non-zero entry in the state)
                    for(typename Vector::iterator it(psi[n][l].begin()); it != psi[n][l].end(); ++it)
                    {
                        // TODO: convert to multiplication operations
                        // Idea: If H0 is diagonal by definition generate a lookup-table or cache.
                        Vector k(h_space_dim);
                        k(get_index_from_iterator(psi[n][l],it)) = value_type(1);
                        value_type Ek = inner_product(k,Ho.apply(g,k));
                        if(Ek != Eo)
                            *it /= (Eo - Ek);
                        else
                            *it = value_type(0);
                    }

                    #ifdef TRACE_STATE_SIZE
                    trace_state_size[n].push_back(psi[n][l].size());
                    #endif //TRACE_STATE_SIZE
                }
                #ifdef EXPANSION_DEBUG
                cout<<"---psi["<<n<<"]["<<l<<"]:---"<<endl;
                psi[n][l].print();
                cout<<"-------------"<<endl;
                #endif //EXPANSION_DEBUG
            }

            #ifdef TRACE_STATE_SIZE
                std::streamsize precision = std::cerr.precision(2);
                std::cerr<<std::endl<<"State Size Trace:"<<std::endl;
                //std::cerr<<"Cluster: Vertices:"<<c.num_vertices()<<"\t Graph:";
                //c.print();
                //std::cerr<<std::endl;
                std::cerr<<"Order \t Avg. size \t % of H-space \t max. size \t % of H-space"<<std::endl;
                for(unsigned int n=1;n<=max_order; ++n)
                {
                    double avg = std::accumulate(trace_state_size[n].begin(),trace_state_size[n].end(),0.)/trace_state_size[n].size();
                    unsigned int max = *std::max_element(trace_state_size[n].begin(),trace_state_size[n].end());
                    double avg_percent = 100*avg / h_space_dim;
                    double max_percent = 100*static_cast<double>(max) / h_space_dim;
                    std::cerr<<n<<"\t"<<static_cast<int>(avg)<<"\t\t";
                    std::cerr<<std::fixed<<avg_percent;
                    std::cerr<<"\t\t"<<max<<"\t\t";
                    std::cerr<<std::fixed<<max_percent<<std::endl;
                }
                std::cerr<<std::endl;
                std::cerr.precision(precision);
            #endif //TRACE_STATE_SIZE
           
            // Convert resulting vector of contribution matrices for different orders
            // to a single matrix with entries that contain all orders -> polynoms
            dense_matrix<Polynomial> result(num_of_degen_states, num_of_degen_states);
            for(int n=max_order; n>=0; --n)
                for(unsigned int j=0; j<num_of_degen_states; ++j)
                    for(unsigned int i=0; i<num_of_degen_states; ++i)
                        result(i,j)[n] = Heff[n](i,j);

            return result;
        }
};

}
#endif //HP2C__ZERO_EXPANSION_HPP
