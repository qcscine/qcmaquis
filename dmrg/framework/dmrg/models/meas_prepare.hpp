/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_MODEL_MEAS_PREPARE_HPP
#define MAQUIS_DMRG_MODEL_MEAS_PREPARE_HPP

#include "dmrg/models/lattice.h"
#include "dmrg/models/generate_mpo.hpp"

namespace meas_prepare {
    
    template<class Matrix, class SymmGroup>
    std::map<std::string, MPO<Matrix,SymmGroup> >
    local(const Lattice & lat,
          std::vector<block_matrix<Matrix, SymmGroup> > const & identities,
          std::vector<block_matrix<Matrix, SymmGroup> > const & fillings,
          std::vector<std::pair<std::vector<block_matrix<Matrix, SymmGroup> >, bool> > const & ops)
	{
        std::map<std::string, MPO<Matrix,SymmGroup> > mpos;
        
        for (std::size_t p = 0; p < lat.size(); ++p)
        {
            if (ops.size() == 1) {
                int type = lat.get_prop<int>("type", p);
                if (ops[0].first[type].n_blocks() > 0) {
                    generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat, identities, fillings);
                    generate_mpo::OperatorTerm<Matrix, SymmGroup> term;
                    term.operators.push_back( std::make_pair(p, ops[0].first[type]) );
                    mpom.add_term(term);
                    
                    mpos[ lat.get_prop<std::string>("label", p) ] = mpom.create_mpo();
                }
            } else {
                std::vector<Lattice::pos_t> neighs = lat.forward(p);
                for (typename std::vector<Lattice::pos_t>::const_iterator hopto = neighs.begin();
                     hopto != neighs.end();
                     ++hopto)
                {
                    int type1 = lat.get_prop<int>("type", p);
                    int type2 = lat.get_prop<int>("type", *hopto);
                    if (ops[0].first[type1].n_blocks() > 0 && ops[1].first[type2].n_blocks() > 0) {
                        generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat, identities, fillings);
                        generate_mpo::OperatorTerm<Matrix, SymmGroup> term;
                        term.operators.push_back( std::make_pair(p, ops[0].first[type1]) );
                        term.operators.push_back( std::make_pair(*hopto, ops[1].first[type2]) );
                        term.with_sign = ops[0].second;
                        mpom.add_term(term);
                        
                        mpos[ lat.get_prop<std::string>("label", p, *hopto) ] = mpom.create_mpo();
                    }
                }
            }
        }
        
        return mpos;
    }
    
    
	template<class Matrix, class SymmGroup>
	MPO<Matrix, SymmGroup>
    average(const Lattice & lat,
            std::vector<block_matrix<Matrix, SymmGroup> > const & identities,
            std::vector<block_matrix<Matrix, SymmGroup> > const & fillings,
            std::vector<std::pair<std::vector<block_matrix<Matrix, SymmGroup> >, bool> > const & ops)
	{
        generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat, identities, fillings);
        
        for (std::size_t p = 0; p < lat.size(); ++p)
        {
            generate_mpo::OperatorTerm<Matrix, SymmGroup> term;
            term.operators.push_back( std::make_pair(p, ops[0].first[lat.get_prop<int>("type", p)]) );
            if (ops.size() == 1) {
                mpom.add_term(term);
            } else {
                term.with_sign = ops[0].second;
            	std::vector<Lattice::pos_t> neighs = lat.forward(p);
            	for (typename std::vector<Lattice::pos_t>::const_iterator hopto = neighs.begin();
            		 hopto != neighs.end();
            		 ++hopto)
            	{
                    generate_mpo::OperatorTerm<Matrix, SymmGroup> term2(term);
                    term2.operators.push_back( std::make_pair(*hopto, ops[1].first[lat.get_prop<int>("type", p)]) );
                    mpom.add_term(term2);
            	}
                
            }
        }
        
        return mpom.create_mpo();
    }
    
}

#endif
