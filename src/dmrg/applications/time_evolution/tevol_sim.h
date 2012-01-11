/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_DMRG_TEVOL_SIM_H
#define APP_DMRG_TEVOL_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/sim/sim.h"

#include "dmrg/sim/te_utils.hpp"
#include "dmrg/mp_tensors/evolve.h"
#include "dmrg/mp_tensors/te.h"

#include "utils/types.h"
//#include "dmrg/utils/noop_storage.h"

using namespace app;

template<class Matrix, class SymmGroup>
struct Trotter_Steps{
	
	/*label type for all DISTINCT unitaries, a list of pairs in format ( [STARTING SITE FOR NON-OVERLAPPING UNITARIES=0 or 1],[TIMESTEP OF UNITARY] )*/
	typedef std::vector< std::pair< std::size_t, typename Matrix::value_type> > u_label_t; 
	
	typedef std::vector<std::map<std::size_t, block_matrix<Matrix, SymmGroup> > > nnop_set_t;
	
	typedef std::vector<typename nnop_set_t::value_type * > diff_seqs_t; // change to boost::ptr_vector
	//typedef diff_seqs_t::iterator diff_seqs_it;
	
	nnop_set_t nnop_set;
	diff_seqs_t diff_seqs;
    
    void clear()
    {
        nnop_set.clear();
        diff_seqs.clear();
    }
};


template <class Matrix, class SymmGroup>
class dmrg_tevol_sim : public sim<Matrix, SymmGroup> {
    
    typedef sim<Matrix, SymmGroup> base;
    
    typedef std::vector<MPOTensor<Matrix, SymmGroup> > mpo_t;
    typedef Boundary<Matrix, SymmGroup> boundary_t;
	
	typedef typename Trotter_Steps<Matrix, SymmGroup>::u_label_t u_label_t;
	typedef typename Matrix::value_type time_step_t; 
	
    enum TEvol_t {te_uknown, te_nn, te_mpo};
	enum Order_t {order_unknown, second_order, fourth_order};
	
public:
    friend
    std::ostream& operator<< (std::ostream& os, TEvol_t const&  t)
    {
        switch (t)
        {
            case te_nn:
                os << "nearest neighbors time-evolve";
                break;
            case te_mpo:
                os << "mpo time-evolve";
                break;
            default:
                os << "uknown time-evolve";
        }
        return os;
    }
	
	friend
    std::ostream& operator<< (std::ostream& os, Order_t const&  o)
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
	
	
    
    dmrg_tevol_sim (DmrgParameters const & parms_, ModelParameters const  & model_)
    : base(parms_, model_, false)
    , parms_orig(parms_)
    , model_orig(model_)
    , te_type( parse_te_type() ) 
    {        
        //trotter_order=parse_order_type();
		
		if (base::parms.template get<std::string>("te_order") == "second")
            trotter_order=second_order;
        else if (base::parms.template get<std::string>("te_order") == "fourth")
            trotter_order=fourth_order;
        else {
            throw std::runtime_error("Don't know this Trotter decomposition");
            trotter_order=order_unknown;
        }
		
		std::cout << "Using " << te_type << std::endl;
		std::cout << "Using " << trotter_order << std::endl;
        
		
		
        base::parms = parms_orig.get_at_index("t", base::sweep);
        base::model = model_orig.get_at_index("t", base::sweep);
        
        base::model_init();
        base::mps_init();
        
        split_H = separate_overlaps(base::H);
        std::cout << split_H.size() << " non overlapping Hamiltonians" << std::endl;
        
        if (te_type == te_nn)
            getUnn(base::parms.template get<double>("dt"),
                         base::sweep < base::parms.template get<int>("nsweeps_img"));
        else if (te_type == te_mpo)
            Umpo = getUmpo(base::parms.template get<double>("dt"),
                           base::sweep < base::parms.template get<int>("nsweeps_img"));
    }
    
    
    int do_sweep (Logger& iteration_log, double time_limit = -1)
    {        
        int pc = 0, mc = 0;
        base::parms = parms_orig.get_at_index("t", base::sweep, &pc);
        base::model = model_orig.get_at_index("t", base::sweep, &mc);
        
        if (mc > 0)
        {
            base::model_init();
            
            split_H = separate_overlaps(base::H);
            std::cout << split_H.size() << " non overlapping Hamiltonians" << std::endl;
            
            if (te_type == te_nn)
                getUnn(base::parms.template get<double>("dt"),
                             base::sweep < base::parms.template get<int>("nsweeps_img"));
            else if (te_type == te_mpo)
                Umpo = getUmpo(base::parms.template get<double>("dt"),
                               base::sweep < base::parms.template get<int>("nsweeps_img"));
        } else {
            if (base::sweep == base::parms.template get<int>("nsweeps_img"))
                if (te_type == te_nn)
                    getUnn(base::parms.template get<double>("dt"), false);
                else if (te_type == te_mpo)
                    Umpo = getUmpo(base::parms.template get<double>("dt"), false);
        }
        
        
        if (te_type == te_nn)
            nn_time_evolve(iteration_log);
        else if (te_type == te_mpo)
            mpo_time_evolve(iteration_log);
        
        double energy = expval(base::mps, base::mpo);
        std::cout << "Energy " << energy << std::endl;
        iteration_log << make_log("Energy", energy);
        
        return -1; // no early exit
    }
    
private:
    
    TEvol_t parse_te_type ()
    {
        if (base::parms.template get<std::string>("te_type") == "nn")
            return te_nn;
        else if (base::parms.template get<std::string>("te_type") == "mpo")
            return te_mpo;
        else {
            throw std::runtime_error("Don't know this type of te");
            return te_uknown;
        }
    }
	Order_t parse_trotter_order ()
    {
        if (base::parms.template get<std::string>("te_order") == "2nd")
            return second_order;
        else if (base::parms.template get<std::string>("te_order") == "4th")
            return fourth_order;
        else {
            throw std::runtime_error("Don't know this Trotter decomposition");
            return order_unknown;
        }
    }
	
	
    /* nearest-neighbors time-evolution */
	
	inline void dist_gates(u_label_t & distinct_gates, time_step_t alpha)
	{
		time_step_t alpha_1,alpha_3;
		
		switch (trotter_order){
			case second_order:
				
				distinct_gates.push_back(std::make_pair (0,alpha*0.5));
				distinct_gates.push_back(std::make_pair (1,alpha));
				
				break;		
			case fourth_order:
				alpha_1=1./(4.0-pow(4.0,0.33333))*alpha;
				alpha_3=alpha-4.0*alpha_1;
				
				distinct_gates.push_back(std::make_pair (0,alpha_1*0.5));
				distinct_gates.push_back(std::make_pair (1,alpha_1));
				distinct_gates.push_back(std::make_pair (0,alpha_3*0.5));
				distinct_gates.push_back(std::make_pair (1,alpha_3));
				
				break;
            default:
                throw std::runtime_error("uknown Trotter decomposition");
                break;

		}
	}
	
	inline void make_seqs(Trotter_Steps<Matrix, SymmGroup> & expH)
	{	
		
		switch (trotter_order){
			case second_order:
				expH.diff_seqs.push_back(&expH.nnop_set[0]);
				expH.diff_seqs.push_back(&expH.nnop_set[1]);
				expH.diff_seqs.push_back(&expH.nnop_set[0]);
				break;
				
			case fourth_order:
				expH.diff_seqs.push_back(&expH.nnop_set[0]);
				expH.diff_seqs.push_back(&expH.nnop_set[1]);
				expH.diff_seqs.push_back(&expH.nnop_set[0]);
				expH.diff_seqs.push_back(&expH.nnop_set[0]);
				expH.diff_seqs.push_back(&expH.nnop_set[1]);
				expH.diff_seqs.push_back(&expH.nnop_set[0]);
				expH.diff_seqs.push_back(&expH.nnop_set[2]);
				expH.diff_seqs.push_back(&expH.nnop_set[3]);
				expH.diff_seqs.push_back(&expH.nnop_set[2]);
				expH.diff_seqs.push_back(&expH.nnop_set[0]);
				expH.diff_seqs.push_back(&expH.nnop_set[1]);
				expH.diff_seqs.push_back(&expH.nnop_set[0]);
				expH.diff_seqs.push_back(&expH.nnop_set[0]);
				expH.diff_seqs.push_back(&expH.nnop_set[1]);
				expH.diff_seqs.push_back(&expH.nnop_set[0]);
				
				break;
            default:
                throw std::runtime_error("uknown Trotter decomposition");
                break;
		}
		
		
	}
	
    
    void
    getUnn (double dt, bool img)
    {
        typename Matrix::value_type I;
        if (img)
            I = utils::real_identity<typename Matrix::value_type>::value;
        else
            I = utils::imag_identity<typename Matrix::value_type>::value;
        typename Matrix::value_type alpha = -I*dt;
        
		/*********/
		
		
		u_label_t distinct_gates;
		
		dist_gates(distinct_gates,alpha);
		
        Unn.clear();
		
		for (typename u_label_t::const_iterator it=distinct_gates.begin(); it!=distinct_gates.end(); ++it){
			
			Unn.nnop_set.push_back(make_exp_nn(split_H[it->first], it->second));
		}
		
		make_seqs(Unn);
				
    }
    
    void nn_time_evolve(Logger& iteration_log)
    {
		
		for(typename Trotter_Steps<Matrix,SymmGroup>::diff_seqs_t::iterator seq_it=Unn.diff_seqs.begin();seq_it!=Unn.diff_seqs.end();seq_it++){
			evolve(base::mps, *(*seq_it), base::parms.template get<std::size_t>("max_bond_dimension"), base::parms.template get<double>("truncation_final"));
		}
	}
    
    
    /* mpo time-evolution */
    
    std::vector<MPO<Matrix, SymmGroup> >
    getUmpo(double dt, bool img)
    {
        typename Matrix::value_type I;
        if (img)
            I = utils::real_identity<typename Matrix::value_type>::value;
        else
            I = utils::imag_identity<typename Matrix::value_type>::value;
        typename Matrix::value_type alpha = -I*dt;
        
        std::vector<MPO<Matrix, SymmGroup> > expMPO(split_H.size(), MPO<Matrix, SymmGroup>(base::lat->size()));
        for (int i=0; i<split_H.size(); ++i)
            expMPO[i] = make_exp_mpo(base::lat->size(), split_H[i], alpha);
        return expMPO;
	}
    
    
    void mpo_time_evolve(Logger& iteration_log)
    {
        for (int which = 0; which < Umpo.size(); ++which)
        {
            time_evolve<Matrix, SymmGroup, NoopStorageMaster> evolution(base::mps,
                                                                        Umpo[which],
                                                                        base::parms, nossm);
            for (int k = 0; k < 5; ++k)
                evolution.sweep(base::sweep, iteration_log);
            evolution.finalize();
            base::mps = evolution.get_current_mps();
        }
        
    }
    
    BaseParameters parms_orig;
    BaseParameters model_orig;
    TEvol_t te_type;
	Order_t trotter_order;
    std::vector<Hamiltonian<Matrix, SymmGroup> > split_H;
    Trotter_Steps<Matrix, SymmGroup> Unn;
    std::vector<MPO<Matrix, SymmGroup> > Umpo;
    
    NoopStorageMaster nossm;
};




#endif
