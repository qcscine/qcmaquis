/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_MODEL_MEAS_PREPARE_HPP
#define MAQUIS_DMRG_MODEL_MEAS_PREPARE_HPP

#include "dmrg/models/lattice.h"
#include "dmrg/models/generate_mpo.hpp"

namespace meas_prepare {
    
    template<class Matrix, class SymmGroup>
    std::pair<std::vector<MPO<Matrix,SymmGroup> >, std::vector<std::string> >
    local(const Lattice & lat,
          block_matrix<Matrix, SymmGroup> const & identity,
          block_matrix<Matrix, SymmGroup> const & fill,
          std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops)
	{
        std::vector<std::string> labels;
        std::vector<MPO<Matrix,SymmGroup> > mpos;
        
        for (std::size_t p = 0; p < lat.size(); ++p)
        {
            if (ops.size() == 1) {
                generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat.size(), identity);
                generate_mpo::Operator_Term<Matrix, SymmGroup> term;
                term.operators.push_back( std::make_pair(p, ops[0].first) );
                term.fill_operator = identity;
                mpom.add_term(term);
                
                mpos.  push_back( mpom.create_mpo()                     );
                labels.push_back( lat.get_prop<std::string>("label", p) );

            } else {
                std::vector<Lattice::pos_t> neighs = lat.forward(p);
                for (typename std::vector<Lattice::pos_t>::const_iterator hopto = neighs.begin();
                     hopto != neighs.end();
                     ++hopto)
                {
                    generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat.size(), identity);
                    generate_mpo::Operator_Term<Matrix, SymmGroup> term;
                    term.operators.push_back( std::make_pair(p, ops[0].first) );
                    term.operators.push_back( std::make_pair(*hopto, ops[1].first) );
                    term.fill_operator = (ops[0].second) ? fill : identity;
                    mpom.add_term(term);

                    mpos.  push_back( mpom.create_mpo()                             );
                    labels.push_back( lat.get_prop<std::string>("label", p, *hopto) );

                }
            }
        }
        
        return std::make_pair( mpos, labels );
    }
    
    
	template<class Matrix, class SymmGroup>
	MPO<Matrix, SymmGroup>
    average(const Lattice & lat,
            block_matrix<Matrix, SymmGroup> const & identity,
            block_matrix<Matrix, SymmGroup> const & fill,
            std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops)
	{
        generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat.size(), identity);
        
        for (std::size_t p = 0; p < lat.size(); ++p)
        {
            generate_mpo::Operator_Term<Matrix, SymmGroup> term;
            term.operators.push_back( std::make_pair(p, ops[0].first) );
            if (ops.size() == 1) {
                mpom.add_term(term);
            } else {
				term.fill_operator = (ops[0].second) ? fill : identity;
            	std::vector<Lattice::pos_t> neighs = lat.forward(p);
            	for (typename std::vector<Lattice::pos_t>::const_iterator hopto = neighs.begin();
            		 hopto != neighs.end();
            		 ++hopto)
            	{
                    generate_mpo::Operator_Term<Matrix, SymmGroup> term2(term);
                    term2.operators.push_back( std::make_pair(*hopto, ops[1].first) );
            	}
                
            }
        }
        
        return mpom.create_mpo();
    }
    
}

#endif
