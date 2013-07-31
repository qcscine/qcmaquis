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

#include "dmrg/sim/te_utils.hpp"
#include "dmrg/mp_tensors/evolve.h"
#include "dmrg/utils/results_collector.h"


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

enum tevol_order_tag {order_unknown, second_order, fourth_order};
inline std::ostream& operator<< (std::ostream& os, tevol_order_tag o)
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


// ******   SIMULATION CLASS   ******
template <class Matrix, class SymmGroup>
class nearest_neighbors_evolver {
public:
    nearest_neighbors_evolver(DmrgParameters * parms_, MPS<Matrix, SymmGroup> * mps_,
                              Lattice_ptr lat_, Hamiltonian<Matrix, SymmGroup> const* H_,
                              int init_sweep=0)
    : parms(parms_)
    , mps(mps_)
    , lat(lat_)
    , H(H_)
    , sweep_(init_sweep)
    , trotter_order(parse_trotter_order((*parms)["te_order"]))
    , block_terms(hamil_to_blocks(*H, lat->size()))
    {
        maquis::cout << "Using nearest-neighbors time evolution." << std::endl;
        maquis::cout << "Using " << trotter_order << std::endl;
        maquis::cout << "Time evolution optimization is "
                     << (((*parms)["te_optim"]) ? "enabled" : "disabled")
                     << std::endl;
        
        /// alpha coeffiecients and set sequence of Uterms according to trotter order
        switch (trotter_order){
			case second_order:
            {
				gates_coeff.push_back(std::make_pair(0,0.5));         // 0-term
				gates_coeff.push_back(std::make_pair(1,1.));          // 1-term
                
                Useq.push_back(0); // odd
				Useq.push_back(1); // even
				Useq.push_back(0); // odd

                if ((*parms)["te_optim"])
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
                
                maquis::cout << "Sequence initialized with " << Useq.size() << " terms." << std::endl;
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
				
                if ((*parms)["te_optim"])
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
                maquis::cout << "Sequence initialized with " << Useq.size() << " terms." << std::endl;

				break;
            }
            default:
            {
                throw std::runtime_error("uknown Trotter decomposition");
                break;
            }
		}
        
        /// compute the time evolution gates
        prepare_te_terms();
    }
    

    void prepare_te_terms()
    {
        size_t L = lat->size();
        
        double dt = (*parms)["dt"];
        typename Matrix::value_type I;
        if (sweep_ < (*parms)["nsweeps_img"])
            I = maquis::traits::real_identity<typename Matrix::value_type>::value;
        else
            I = maquis::traits::imag_identity<typename Matrix::value_type>::value;
        typename Matrix::value_type alpha = -I * dt;
        
        Uterms.resize(gates_coeff.size(), trotter_gate<Matrix, SymmGroup>(L));
        for (size_t i=0; i<gates_coeff.size(); ++i) {
            Uterms[i].clear();
            Uterms[i].pfirst = gates_coeff[i].first;
            for (size_t p=gates_coeff[i].first; p<L-1; p+=2){
                if ((*parms)["expm_method"] == "heev")
                    Uterms[i].add_term(p, op_exp_hermitian(H->get_phys()*H->get_phys(), block_terms[p], gates_coeff[i].second*alpha));
                else
                    Uterms[i].add_term(p, op_exp(H->get_phys()*H->get_phys(), block_terms[p], gates_coeff[i].second*alpha));
            }
        }
    }
    
    void operator()(int nsteps)
    {
        iteration_results_.clear();
        
        int ns = sweep_ + nsteps;
        
        if (nsteps < 2 || !(*parms)["te_optim"]) {
            // nsteps sweeps
            for (int i=sweep_; i < ns; ++i)
            {
                sweep_ = i;
                (*parms).set("sweep", sweep_);
                evolve_time_step(Useq);
            }
        } else {
            // one sweep
            (*parms).set("sweep", sweep_);
            evolve_time_step(Useq_bmeas);
            
            // nsteps - 2 sweeps
            for (int i=sweep_+1; i<ns-1; ++i) {
                sweep_ = i;
                (*parms).set("sweep", sweep_);
                evolve_time_step(Useq_double);
            }
            
            // one sweep
            sweep_ += 1;
            (*parms).set("sweep", sweep_);
            evolve_time_step(Useq_ameas);
        }
        assert(sweep_ == ns-1);
    }
    
    int sweep() const
    {
        return sweep_;
    }
    
    results_collector const& iteration_results() const
    {
        return iteration_results_;
    }
    
private:
    tevol_order_tag parse_trotter_order (std::string const & trotter_param)
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

    void evolve_time_step(std::vector<std::size_t> const & gates_i)
    {
        assert(gates_i.size() > 0);
        
        for (size_t i=0; i<gates_i.size(); ++i) {
            if (mps->canonization(true) < mps->length()/2)
                evolve_l2r(gates_i[i]);
            else
                evolve_r2l(gates_i[i]);
        }
    }
    
    void evolve_l2r(std::size_t gate_index)
    {
        std::size_t L = mps->length();
        std::vector<block_matrix<Matrix, SymmGroup> > const & ops = Uterms[gate_index].vgates;
        std::vector<long> const & idx = Uterms[gate_index].idx;
        int pfirst = Uterms[gate_index].pfirst;
        std::size_t Mmax=(*parms)["max_bond_dimension"];
        double cutoff=(*parms)["truncation_final"];
        MPS<Matrix, SymmGroup> const& constmps = *mps;
        
        assert(L == idx.size());

        if (mps->canonization() != pfirst + 1)
            mps->canonize(pfirst + 1);
        
        // TODO: remove pfirst!
        for (std::size_t p = pfirst; p <= L-1; p += 2)
        {
            if (idx[p] != -1)
            {
                constmps[p].make_left_paired();
                constmps[p+1].make_right_paired();
                
                block_matrix<Matrix, SymmGroup> v0, v1;
                gemm(constmps[p].data(), constmps[p+1].data(), v0); // outer product of two sites
                
                v1 = contraction::multiply_with_twosite(v0, ops[idx[p]],
                                                        constmps[p].row_dim(), constmps[p+1].col_dim(),
                                                        constmps[p].site_dim());
                truncation_results trunc = compression::replace_two_sites_l2r(*mps, Mmax, cutoff, v1, p);
                iteration_results_["BondDimension"]     << trunc.bond_dimension;
                iteration_results_["TruncatedWeight"]   << trunc.truncated_weight;
                iteration_results_["TruncatedFraction"] << trunc.truncated_fraction;
                iteration_results_["SmallestEV"]        << trunc.smallest_ev;
            }
            mps->move_normalization_l2r(p+1, p+3, DefaultSolver());
        }
        mps->canonization(true);
        assert(mps->canonization() == L-1);
        // maquis::cout << "Norm loss " << i << ": " << trace(t) << " " << -log(trace(t)) << std::endl;
    }
    
    void evolve_r2l(std::size_t gate_index)
    {
        std::size_t L = mps->length();
        std::vector<block_matrix<Matrix, SymmGroup> > const & ops = Uterms[gate_index].vgates;
        std::vector<long> const & idx = Uterms[gate_index].idx;
        int pfirst = Uterms[gate_index].pfirst;
        std::size_t Mmax=(*parms)["max_bond_dimension"];
        double cutoff=(*parms)["truncation_final"];
        MPS<Matrix, SymmGroup> const& constmps = *mps;

        assert(L == idx.size());
        
        int startpos = std::min(L-2-(L-pfirst)%2, L-2);
        if (mps->canonization() != startpos)
            mps->canonize(startpos);
        
        for (int p = std::min(L-2-(L-pfirst)%2, L-2); p >= pfirst; p -= 2)
        {
            if (idx[p] != -1)
            {
                constmps[p].make_left_paired();
                constmps[p+1].make_right_paired();
                
                block_matrix<Matrix, SymmGroup> v0, v1;
                gemm(constmps[p].data(), constmps[p+1].data(), v0); // outer product of two sites
                
                v1 = contraction::multiply_with_twosite(v0, ops[idx[p]],
                                                        constmps[p].row_dim(), constmps[p+1].col_dim(),
                                                        constmps[p].site_dim());
                truncation_results trunc = compression::replace_two_sites_r2l(*mps, Mmax, cutoff, v1, p);
                iteration_results_["BondDimension"]     << trunc.bond_dimension;
                iteration_results_["TruncatedWeight"]   << trunc.truncated_weight;
                iteration_results_["TruncatedFraction"] << trunc.truncated_fraction;
                iteration_results_["SmallestEV"]        << trunc.smallest_ev;
            }
            mps->move_normalization_r2l(p, std::max(static_cast<long>(p)-2,0L), DefaultSolver());
        }
        
        
        mps->canonization(true);
        assert(mps->canonization() == 0);
        // maquis::cout << "Norm loss " << i << ": " << trace(t) << " " << -log(trace(t)) << std::endl;
    }

private:
    DmrgParameters * parms;
    MPS<Matrix, SymmGroup> * mps;
    Lattice_ptr lat;
    Hamiltonian<Matrix, SymmGroup> const * H;
    int sweep_;
    
    results_collector iteration_results_;
    
	tevol_order_tag trotter_order;
    std::vector<block_matrix<Matrix, SymmGroup> > block_terms;
    std::vector<std::pair<std::size_t,double> > gates_coeff;
    std::vector<trotter_gate<Matrix, SymmGroup> > Uterms;
    std::vector<std::size_t> Useq; // trivial sequence
    std::vector<std::size_t> Useq_double, Useq_bmeas, Useq_ameas; // sequence with two sweeps; meas before; meas after
};

#endif
