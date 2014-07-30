/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MODELS_CODED_NONE_H
#define MODELS_CODED_NONE_H

#include <sstream>

#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"

/* ****************** BOSE-HUBBARD */
template<class Matrix>
class BoseHubbardNone : public model_impl<Matrix, TrivialGroup>
{
    typedef model_impl<Matrix, TrivialGroup> base;

    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    
    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;
    
    typedef typename base::size_t size_t;
    typedef typename Matrix::value_type value_type;
    
public:
    BoseHubbardNone (const Lattice& lat, BaseParameters & model_)
    : model(model_)
    , lattice(lat)
    , tag_handler(new table_type())
    {
        int Nmax = model["Nmax"];
        double U = model["U"];
        double t = model["t"];
        double V = model["V"];
        
        op_t ident_op;
        op_t create_op, destroy_op, count_op, interaction_op;

        TrivialGroup::charge C = TrivialGroup::IdentityCharge;
        size_t N = Nmax+1;
        
        phys.insert(std::make_pair(C, N));
        ident_op.insert_block(Matrix::identity_matrix(N), C, C);
        
        Matrix mcount(N,N), minteraction(N,N), mcreate(N,N), mdestroy(N,N);
        for (int n=1; n<=Nmax; ++n)
        {
            mcount(n,n) = n;
            if ((n*n-n) != 0)
                minteraction(n,n) = n*n-n;
            
            mcreate(n-1,n) = std::sqrt(value_type(n));   // input n-1, output n
            mdestroy(n,n-1) = std::sqrt(value_type(n));  // input n,   output n-1
        }
        count_op.insert_block(mcount, C,C);
        interaction_op.insert_block(minteraction, C,C);
        create_op.insert_block(mcreate, C,C);
        destroy_op.insert_block(mdestroy, C,C);
        
        /**********************************************************************/
        /*** Create operator tag table ****************************************/
        /**********************************************************************/
        
#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);
        
        REGISTER(ident,   tag_detail::bosonic)
        REGISTER(create,  tag_detail::bosonic)
        REGISTER(destroy, tag_detail::bosonic)
        REGISTER(count,   tag_detail::bosonic)
        REGISTER(interaction,   tag_detail::bosonic)
        
#undef REGISTER
        /**********************************************************************/

        
        for (int p=0; p<lat.size(); ++p) {
            /* interaction */
            if (U != 0.) {
                term_descriptor term;
                term.is_fermionic = false;
                term.coeff = U/2.;
                term.push_back( boost::make_tuple(p, interaction) );
                this->terms_.push_back(term);
            }
            
            std::vector<int> neighs = lat.forward(p);
            for (int n=0; n<neighs.size(); ++n) {
                /* hopping */
                {
                    term_descriptor term;
                    term.coeff = -t;
                    term.push_back( boost::make_tuple(p, create) );
                    term.push_back( boost::make_tuple(neighs[n], destroy) );
                    this->terms_.push_back(term);
                }
                {
                    term_descriptor term;
                    term.coeff = -t;
                    term.push_back( boost::make_tuple(p, destroy) );
                    term.push_back( boost::make_tuple(neighs[n], create) );
                    this->terms_.push_back(term);
                }
                /* nearest-neighborn interaction */
                if (V != 0.){
                    term_descriptor term;
                    term.coeff = V;
                    term.push_back( boost::make_tuple(p, count) );
                    term.push_back( boost::make_tuple(neighs[n], count) );
                    this->terms_.push_back(term);
                }
            }
        }
        
    }
    
    void update(BaseParameters const& p)
    {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
        return;
    }
    
    Index<TrivialGroup> const& phys_dim(size_t type) const
    {
        return phys;
    }
    tag_type identity_matrix_tag(size_t type) const
    {
        return ident;
    }
    tag_type filling_matrix_tag(size_t type) const
    {
        return identity_matrix_tag(type);
    }
    typename TrivialGroup::charge total_quantum_numbers(BaseParameters & parms) const
    {
        return TrivialGroup::IdentityCharge;
    }
    
    
   measurements_type measurements() const
    {
        typedef std::vector<block_matrix<Matrix, TrivialGroup> > op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
        
        measurements_type meas;
        
        if (model["MEASURE[Density]"]) {
            meas.push_back( new measurements::average<Matrix, TrivialGroup>("Density", lattice,
                                                                  op_vec(1,this->identity_matrix(0)),
                                                                  op_vec(1,this->filling_matrix(0)),
                                                                  op_vec(1,tag_handler->get_op(count))) );
        }
        
        if (model["MEASURE[Local density]"]) {
            meas.push_back( new measurements::local<Matrix, TrivialGroup>("Local density", lattice,
                                                                op_vec(1,this->identity_matrix(0)),
                                                                op_vec(1,this->filling_matrix(0)),
                                                                op_vec(1,tag_handler->get_op(count))) );
        }
        
        if (model["MEASURE[Onebody density matrix]"]) {
            bond_element ops;
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(create)), false) );
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(destroy)), false) );
            meas.push_back( new measurements::correlations<Matrix, TrivialGroup>("Onebody density matrix", lattice,
                                                                                 op_vec(1,this->identity_matrix(0)),
                                                                                 op_vec(1,this->filling_matrix(0)),
                                                                                 ops, true, false) );
        }
        
        return meas;
    }
    
    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "n")
            return count;
        else if (name == "bdag")
            return create;
        else if (name == "b")
            return destroy;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return 0;
    }
    
    table_ptr operators_table() const
    {
        return tag_handler;
    }

    
private:
    BaseParameters & model;
    Lattice lattice;
    Index<TrivialGroup> phys;

    table_ptr tag_handler;
    tag_type ident, create, destroy, count, interaction;
};



#endif
