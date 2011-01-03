#ifndef HP2C__HIGH_T_EXPANSION_HPP
#define HP2C__HIGH_T_EXPANSION_HPP

#include "numeric/inner_product.hpp"
#include <cmath>


namespace series_expansion
{

template <typename Operator, typename Graph, typename Vector, typename Polynomial>
class high_t_expansion
{
    private:
        Operator H;
        unsigned int max_order;
        Polynomial expanded_log(Polynomial const& t, unsigned int order) const
        {
            // TODO implement general formula for this
            Polynomial result;
            result[0] = 0;
            result[1] = t[1]/t[0];
            result[2] = (-t[1]*t[1] + 2*t[0]*t[2])/(2*t[0]*t[0]);
            result[3] = ( t[1]*t[1]*t[1] - 3*t[0]*t[1]*t[2] + 3*t[0]*t[0]*t[3])/(3*t[0]*t[0]*t[0]);
            result[4] = ( -t[1]*t[1]*t[1]*t[1] + 4*t[0]*t[1]*t[1]*t[2] - 2*t[0]*t[0]*t[2]*t[2] - 4*t[0]*t[0]*t[1]*t[3] + 4*t[0]*t[0]*t[0]*t[4])/(4*t[0]*t[0]*t[0]*t[0]);
            result[5] = ( t[1]*t[1]*t[1]*t[1]*t[1] - 5*t[0]*t[1]*t[1]*t[1]*t[2] + 5*t[0]*t[0]*t[1]*t[2]*t[2] + 5*t[0]*t[0]*t[1]*t[1]*t[3] - 5*t[0]*t[0]*t[0]*t[2]*t[3] - 5*t[0]*t[0]*t[0]*t[1]*t[4] + 5*t[0]*t[0]*t[0]*t[0]*t[5])/(5*t[0]*t[0]*t[0]*t[0]*t[0]);
            result[6] = ( -t[1]*t[1]*t[1]*t[1]*t[1]*t[1] + 6*t[0]*t[1]*t[1]*t[1]*t[1]*t[2] - 9*t[0]*t[0]*t[1]*t[1]*t[2]*t[2] + 2*t[0]*t[0]*t[0]*t[2]*t[2]*t[2] - 6*t[0]*t[0]*t[1]*t[1]*t[1]*t[3] + 12*t[0]*t[0]*t[0]*t[1]*t[2]*t[3] - 3*t[0]*t[0]*t[0]*t[0]*t[3]*t[3]  + 6*t[0]*t[0]*t[0]*t[1]*t[1]*t[4] - 6*t[0]*t[0]*t[0]*t[0]*t[2]*t[4] - 6*t[0]*t[0]*t[0]*t[0]*t[1]*t[5] + 6*t[0]*t[0]*t[0]*t[0]*t[0]*t[6])/(6*t[0]*t[0]*t[0]*t[0]*t[0]*t[0]);
            result.truncate(order);
            return result;
        }

    public:
        typedef typename Vector::value_type value_type;

        high_t_expansion(Operator const& H, unsigned int max_order)
            : H(H), max_order(max_order)
            {
                BOOST_STATIC_ASSERT(( boost::is_same<typename Vector::value_type, typename Polynomial::value_type>::value ));
            }

        ~high_t_expansion()
            {}
	
        std::pair<Polynomial,std::vector<value_type> > operator() (Graph const& g) const
        {
            typename Vector::size_type h_space_dim = calc_hilbert_space_dimension(typename Vector::size_type(),H,g);

            #ifdef TRACE_STATE_SIZE
                // This flag enables a counter that counts how many entries the sparse vector has
                // and prints out the percentage of the full hilbert space dimension.
                std::vector<std::vector<unsigned int> > trace_state_size(max_order/2+1);
            #endif //TRACE_STATE_SIZE

            // Compute trace of the Operator Ho.
            // ( Tr(Ho) )
            Polynomial result;
            for(typename Vector::size_type i(0); i != h_space_dim; ++i)
            {
                Vector state(h_space_dim,value_type(0));
                state(i) = value_type(1);
                
                // 0th order
                assert( inner_product(state,state) == value_type(1) );
                result[0] += value_type(1);
                
                // 1st - <order>th order
                //
                // we can half of the computations, since we can
                // evolve the state from both sides:
                // result[4] = inner_product(<init|*Ho*Ho ,Ho* Ho*|init> )
                for(unsigned int n=1; n<= max_order/2; ++n)
                {
                    Vector state2(state);
                    state = H.apply(g,state);
                    #ifdef TRACE_STATE_SIZE
                        trace_state_size[n].push_back(state.size());
                    #endif //TRACE_STATE_SIZE
                    result[2*n] += inner_product(state,state);
                    result[2*n-1] += inner_product(state,state2);
                }
            }

            #ifdef TRACE_STATE_SIZE
            std::streamsize precision = std::cerr.precision(2);
            std::cerr<<std::endl<<"State Size Trace:"<<std::endl;
            //std::cerr<<"Cluster: Vertices:"<<c.num_vertices()<<"\t Graph:";
            //c.print();
            //std::cerr<<std::endl;
            std::cerr<<"Order \t Avg. size \t % of H-space \t max. size \t % of H-space"<<std::endl;
            
            for(unsigned int n=1;n<=max_order/2; ++n)
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

            // TODO convert that to multiplications by integers
            // Devide each term by 1/(n!) , where n is the order
            value_type fac(1);
            std::vector<value_type> denominators(max_order+1);
            denominators[0] = 1;
            for(unsigned int i=1; i<=max_order; ++i)
            {
                fac*=i;
                denominators[i] = fac * std::pow(static_cast<unsigned int>(Operator::normalization_factor),i);
            }

            //result = expanded_log(result,max_order); // t[0] ist Anzahl Zustände = d^N, nachdenken ....
            return make_pair(result,denominators);
        }

};

}
#endif //HP2C__HIGH_T_EXPANSION_HPP
