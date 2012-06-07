/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_DMRG_TEVOL_NN_SIM_H
#define APP_DMRG_TEVOL_NN_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include "tevol_sim.h"

#include "dmrg/sim/te_utils.hpp"
#include "dmrg/mp_tensors/evolve.h"

#include "utils/types.h"


// ******     HELPER OBJECTS     ******
template<class Matrix, class SymmGroup>
struct trotter_gate {
    typedef std::vector<long> idx_t;
    typedef std::vector<block_matrix<Matrix, SymmGroup> > vgates_t;
    
    std::size_t pfirst;
	idx_t idx;
	vgates_t vgates;
    
    trotter_gate(std::size_t L) : idx(L, -1) { }
    
    void add_term(std::size_t p, block_matrix<Matrix, SymmGroup> const & block)
    {
        assert(idx[p] == -1);
        vgates.push_back(block);
        idx[p] = vgates.size()-1;
    }
    
    void clear()
    {
        vgates.clear();
        idx = idx_t(idx.size(), -1);
        pfirst = 0;
    }
};


// ******   SIMULATION CLASS   ******
template <class Matrix, class SymmGroup>
class dmrg_tevol_nn_sim : public dmrg_tevol_sim<Matrix, SymmGroup> {
    
    typedef dmrg_tevol_sim<Matrix, SymmGroup> base;
    	
	enum order_t {order_unknown, second_order, fourth_order};
	
public:
    typedef typename base::mpo_t mpo_t;
    typedef typename base::boundary_t boundary_t;

	friend
    std::ostream& operator<< (std::ostream& os, order_t const&  o)
    {
        switch (o)
        {
            case second_order:
                os << "Second order Trotter decomposition";
                break;
            case fourth_order:
                os << "Fourth order Trotter decomposition";
                break;
            default:
                os << "uknown Trotter decomposition";
        }
        return os;
    }
	
    dmrg_tevol_nn_sim (DmrgParameters const & parms_, ModelParameters const  & model_)
    : base(parms_, model_)
    , trotter_order(parse_trotter_order(this->parms.template get<std::string>("te_order")))
    {
        maquis::cout << "Using nearest-neighbors time evolution." << std::endl;
        maquis::cout << "Using " << trotter_order << std::endl;
        maquis::cout << "Time evolution optimization is "
                     << ((this->parms.template get<bool>("te_optim")) ? "enabled" : "disabled")
                     << std::endl;
        
        // alpha coeffiecients and set sequence of Uterms according to trotter order
        switch (trotter_order){
			case second_order:
            {
				gates_coeff.push_back(std::make_pair(0,0.5));         // 0-term
				gates_coeff.push_back(std::make_pair(1,1.));          // 1-term
                
                Useq.push_back(0); // odd
				Useq.push_back(1); // even
				Useq.push_back(0); // odd

                if (this->parms.template get<bool>("te_optim"))
                {
                    gates_coeff.push_back(std::make_pair(0,1.));     // 2-term

                    Useq_bmeas.push_back(0); // odd
                    Useq_bmeas.push_back(1); // even
                    Useq_bmeas.push_back(2); // odd

                    Useq_double.push_back(1); // even
                    Useq_double.push_back(2); // odd
                    
                    Useq_ameas.push_back(1); // even
                    Useq_ameas.push_back(0); // odd
                }
                
                maquis::cout << "Sequence initialiez with " << Useq.size() << " terms." << std::endl;
				break;		
            }
            case fourth_order:
			{
                double alpha_1=1./(4.0-pow(4.0,0.33333));
				double alpha_3=1.-4.0*alpha_1;
				
				gates_coeff.push_back(std::make_pair(0,alpha_1*0.5)); // 0-term
				gates_coeff.push_back(std::make_pair(1,alpha_1));     // 1-term
				gates_coeff.push_back(std::make_pair(0,alpha_3*0.5)); // 2-term
				gates_coeff.push_back(std::make_pair(1,alpha_3));     // 3-term
                
				Useq.push_back(0); // odd
				Useq.push_back(1); // even
				Useq.push_back(0); // odd
				Useq.push_back(0); // odd
				Useq.push_back(1); // even
				Useq.push_back(0); // odd
				Useq.push_back(2); // odd
				Useq.push_back(3); // even
				Useq.push_back(2); // odd
				Useq.push_back(0); // odd
				Useq.push_back(1); // even
				Useq.push_back(0); // odd
				Useq.push_back(0); // odd
				Useq.push_back(1); // even
				Useq.push_back(0); // odd
				
                if (this->parms.template get<bool>("te_optim"))
                {
                    gates_coeff.push_back(std::make_pair(0,alpha_1));                 // 4-term
                    gates_coeff.push_back(std::make_pair(0,alpha_3*0.5+alpha_1*0.5)); // 5-term
                    
                    Useq_bmeas.push_back(0); // odd
                    Useq_bmeas.push_back(1); // even
                    Useq_bmeas.push_back(4); // odd
                    Useq_bmeas.push_back(1); // even
                    Useq_bmeas.push_back(5); // odd
                    Useq_bmeas.push_back(3); // even
                    Useq_bmeas.push_back(5); // odd
                    Useq_bmeas.push_back(1); // even
                    Useq_bmeas.push_back(4); // odd
                    Useq_bmeas.push_back(1); // even
                    Useq_bmeas.push_back(4); // odd
                    
                    Useq_double.push_back(1); // even
                    Useq_double.push_back(4); // odd
                    Useq_double.push_back(1); // even
                    Useq_double.push_back(5); // odd
                    Useq_double.push_back(3); // even
                    Useq_double.push_back(5); // odd
                    Useq_double.push_back(1); // even
                    Useq_double.push_back(4); // odd
                    Useq_double.push_back(1); // even
                    Useq_double.push_back(4); // odd
                    
                    Useq_ameas.push_back(1); // even
                    Useq_ameas.push_back(4); // odd
                    Useq_ameas.push_back(1); // even
                    Useq_ameas.push_back(5); // odd
                    Useq_ameas.push_back(3); // even
                    Useq_ameas.push_back(5); // odd
                    Useq_ameas.push_back(1); // even
                    Useq_ameas.push_back(4); // odd
                    Useq_ameas.push_back(1); // even
                    Useq_ameas.push_back(0); // odd
                }
                maquis::cout << "Sequence initialied with " << Useq.size() << " terms." << std::endl;

				break;
            }
            default:
            {
                throw std::runtime_error("uknown Trotter decomposition");
                break;
            }
		}
    }
    

protected:
    void prepare_te_terms(bool split_hamil=true)
    {
        size_t L = this->lat->size();
        if (split_hamil)
            block_terms = hamil_to_blocks(this->H, L);
        
        typename Matrix::value_type I;
        if (this->sweep < this->parms.template get<int>("nsweeps_img"))
            I = utils::real_identity<typename Matrix::value_type>::value;
        else
            I = utils::imag_identity<typename Matrix::value_type>::value;
        typename Matrix::value_type alpha = -I*this->parms.template get<double>("dt");
        
        Uterms.resize(gates_coeff.size(), trotter_gate<Matrix, SymmGroup>(L));
        for (size_t i=0; i<gates_coeff.size(); ++i) {
            Uterms[i].clear();
            Uterms[i].pfirst = gates_coeff[i].first;
            for (size_t p=gates_coeff[i].first; p<L-1; p+=2){
                Uterms[i].add_term(p, op_exp(this->H.get_phys()*this->H.get_phys(), block_terms[p], gates_coeff[i].second*alpha));
            }
        }
    }
    
    void evolve_time_step(Logger & iteration_log)
    { return evolve_time_step(Useq, iteration_log); }
    
    void evolve_time_step(std::vector<std::size_t> const & gates_i, Logger & iteration_log)
    {
        assert(gates_i.size() > 0);
        for (size_t i=0; i<gates_i.size(); ++i) {
            if (this->mps.canonization(true) < this->mps.length()/2)
                evolve_l2r(this->mps, Uterms[gates_i[i]].vgates, Uterms[gates_i[i]].idx, Uterms[gates_i[i]].pfirst,
                           this->parms.template get<std::size_t>("max_bond_dimension"),
                           this->parms.template get<double>("truncation_final"),
                           &iteration_log);
            else
                evolve_r2l(this->mps, Uterms[gates_i[i]].vgates, Uterms[gates_i[i]].idx, Uterms[gates_i[i]].pfirst,
                           this->parms.template get<std::size_t>("max_bond_dimension"),
                           this->parms.template get<double>("truncation_final"),
                           &iteration_log);
        }
    }
    
    void evolve_ntime_steps(int nsteps, Logger & iteration_log)
    {
        int ns = this->sweep + nsteps;

        if (nsteps < 2 || !this->parms.template get<bool>("te_optim")) {
            // nsteps sweeps
            for (; this->sweep < ns; ++(this->sweep))
            {
                this->parms.set("sweep", this->sweep);
                evolve_time_step(Useq, iteration_log);
            }
        } else {
            this->parms.set("sweep", this->sweep);
            // one sweep
            evolve_time_step(Useq_bmeas, iteration_log);
            ++(this->sweep);
            
            // nsteps - 2 sweep
            for (; this->sweep < ns-1; ++(this->sweep))
                evolve_time_step(Useq_double, iteration_log);
            
            // one sweep
            evolve_time_step(Useq_ameas, iteration_log);
            ++(this->sweep);
        }
        // this->sweep = ns!
        --(this->sweep);
    }

    
private:
    
	order_t parse_trotter_order (std::string const & trotter_param)
    {
        if (trotter_param == "second")
            return second_order;
        else if (trotter_param == "fourth")
            return fourth_order;
        else {
            throw std::runtime_error("Don't know this Trotter decomposition");
            return order_unknown;
        }
    }
    
        
	order_t trotter_order;
    std::vector<block_matrix<Matrix, SymmGroup> > block_terms;
    std::vector<std::pair<std::size_t,double> > gates_coeff;
    std::vector<trotter_gate<Matrix, SymmGroup> > Uterms;
    std::vector<std::size_t> Useq; // trivial sequence
    std::vector<std::size_t> Useq_double, Useq_bmeas, Useq_ameas; // sequence with two sweeps; meas before; meas after
};

#endif
