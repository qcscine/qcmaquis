#ifndef HP2C__OPERATORS_HPP
#define HP2C__OPERATORS_HPP

//#include <state.hpp>
#include "numeric/get_index_from_iterator.hpp"
#include "numeric/index_from_site_state.hpp"
#include <utility>

namespace series_expansion
{

//
// unperturbed (diagonal) part of the Hamiltonian H = H0 + lambda * V
//

class heisenberg_h0
{
	public:
        enum {singlet, trp_minus, trp_null, trp_plus };
        enum {states_per_site = trp_plus + 1 }; // StatesPerSite
        enum {SPS = states_per_site };          // just for convinience (see above)
        enum {normalization_factor = 4};        // in principle we are want 1/4 and 3/4, this is the denominator

	    template <typename Cluster, typename Vector>	
        const Vector apply(const Cluster& c,const Vector &s) const
		{	
            typedef typename Vector::value_type value_type;
			Vector result(s);
			for(typename Vector::iterator it(result.begin()); it != result.end(); ++it)
			{
                if( *it == value_type(0) ) continue; // We are only interested in non-zero entries,
                                                    // TODO: Perhaps we should write a
                                                    // return_next_non-zero_entry() function
                
                // TODO this is not generic (arbitrary vector type)
				// Count the SINGELTs and TRIPLETs for each basis state and multiply the amplitude with the result
                value_type a(0);
				for(unsigned int i=0; i<num_vertices(c); ++i)
				{
					const unsigned int sstate = 
                        get_site_state_from_index<SPS>(
                                get_index_from_iterator(result,it),
                                i);
					if(sstate == singlet)
					{
						a -= value_type(3);
					}
					else
					{
						a += value_type(1);
					}
				}
				*it *= a;
			}
			return result;
		}
};


//
// Perturbation of the Hamiltonian H= H0 + lambda * V
//
class heisenberg_v
{
	public:
        enum {singlet, trp_minus, trp_null, trp_plus };
        enum {states_per_site = trp_plus + 1 }; // StatesPerSite
        enum {SPS = states_per_site };           // just for convinience (see above)
        enum {normalization_factor = 4};        // in principle we are want 1/4 and 3/4, this is the denominator

        template <typename Cluster, typename Vector>
		const Vector apply (const Cluster& c, const Vector &s) const
		{
            typedef typename Vector::value_type value_type;
			Vector result(s.size());
			for(typename Vector::const_iterator it(s.begin()); it != s.end(); ++it)
			{
                if( *it == value_type(0) ) continue; // We are only interested in non-zero entries,
                                                    // TODO: Perhaps we should write a
                                                    // return_next_non-zero_entry() function
                
                const typename Vector::size_type i = get_index_from_iterator(s,it);
				const value_type cf(*it * value_type(2));
				
                typename Cluster::edge_iterator c_it, c_end;
				for(boost::tie(c_it, c_end) = edges(c); c_it != c_end; ++c_it)
				{
					const unsigned int siteA = source(*c_it,c);
					const unsigned int siteB = target(*c_it,c);
					const unsigned int stateA = get_site_state_from_index<SPS>(i,siteA);
					const unsigned int stateB = get_site_state_from_index<SPS>(i,siteB);
					
                    //
                    // Operations checked on 2010 12 15
                    // 
                    if(stateA == singlet && stateB == singlet)
                    {
                        result(index_from_site_state_pair<SPS>(i,siteA,trp_plus,siteB,trp_minus))   += -cf;
                        result(index_from_site_state_pair<SPS>(i,siteA,trp_minus,siteB,trp_plus))   += -cf;
                        result(index_from_site_state_pair<SPS>(i,siteA,trp_null,siteB,trp_null))    += cf;

                    }
                    else if ( stateA == singlet || stateB == singlet) // XOR (due to previous if (.. &&..) )
                    {
                        result(index_from_site_state_pair<SPS>(i,siteA,stateB,siteB,stateA))       += cf;
                    }
                    else if ( (stateA == trp_plus && stateB == trp_plus) || (stateA == trp_minus && stateB == trp_minus) )
                    {
                        result(i)   += cf;
                    }
                    else if ( (stateA == trp_plus && stateB == trp_minus) || (stateA == trp_minus && stateB == trp_plus) )
                    {
                        result(i)   += -cf;

                        result(index_from_site_state_pair<SPS>(i,siteA,trp_null,siteB,trp_null))    += cf;
                       
                        result(index_from_site_state_pair<SPS>(i,siteA,singlet,siteB,singlet))      += -cf;
                    }
                    else if ( stateA == trp_null && stateB == trp_null)
                    {
                        result(index_from_site_state_pair<SPS>(i,siteA,trp_plus,siteB,trp_minus))   += cf;
                       
                        result(index_from_site_state_pair<SPS>(i,siteA,trp_minus,siteB,trp_plus))   += cf;

                        result(index_from_site_state_pair<SPS>(i,siteA,singlet,siteB,singlet))      += cf;
                    }
                    else if (stateA == trp_null || stateB == trp_null ) // XOR due to previous if( .. && .. )
                    {
                        result(index_from_site_state_pair<SPS>(i,siteA,stateB,siteB,stateA))        += cf;
                    }
				}
			}
			return result;
		}
};

class heisenberg_hamiltonian
{
	public:
        enum {singlet, trp_minus, trp_null, trp_plus };
        enum {states_per_site = trp_plus + 1 }; // StatesPerSite
        enum {normalization_factor = 4};        // in principle we are want 1/4 and 3/4, this is the denominator

        heisenberg_h0 Ho;
        heisenberg_v  V;

        template <typename Cluster, typename Vector>
		const Vector apply (const Cluster& c, const Vector &s) const
		{
            Vector result(Ho.apply(c,s));
            result += (V.apply(c,s));
            return result;
        }
};

}

#endif //HP2C__OPERATORS_HPP
