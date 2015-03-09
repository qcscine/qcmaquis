/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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

#ifndef DMRG_CONTINUOUS_MODELS_U1_H
#define DMRG_CONTINUOUS_MODELS_U1_H

#include <sstream>

#include <cmath>

#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "dmrg/utils/BaseParameters.h"

/* ****************** OPTICAL LATTICE (with symmetry) */
template<class Matrix>
class OpticalLattice : public model_impl<Matrix, U1>
{
    typedef model_impl<Matrix, U1> base;
    
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    
    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;
    
    typedef typename Matrix::value_type value_type;
public:
    OpticalLattice (const Lattice& lat_, BaseParameters & model_)
    : lat(lat_)
    , model(model_)
    , tag_handler(new table_type())
    {
        int Nmax = model["Nmax"];
        
        op_t ident_op;
        op_t create_op, destroy_op, count_op, interaction_op;
        
        phys.insert(std::make_pair(0, 1));
        ident_op.insert_block(Matrix(1, 1, 1), 0, 0);
        
        for (int n=1; n<=Nmax; ++n)
        {
            phys.insert(std::make_pair(n, 1));
            
            ident_op.insert_block(Matrix(1, 1, 1), n, n);
            
            count_op.insert_block(Matrix(1, 1, n), n, n);
            if ((n*n-n) != 0)
                interaction_op.insert_block(Matrix(1, 1, n*n-n), n, n);
            
            
            create_op.insert_block(Matrix(1, 1, std::sqrt(value_type(n))), n-1, n);
            destroy_op.insert_block(Matrix(1, 1, std::sqrt(value_type(n))), n, n-1);
        }
        
        
        /**********************************************************************/
        /*** Create operator tag table ****************************************/
        /**********************************************************************/
        
#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);
        
        REGISTER(ident,       tag_detail::bosonic)
        REGISTER(create,      tag_detail::bosonic)
        REGISTER(destroy,     tag_detail::bosonic)
        REGISTER(count,       tag_detail::bosonic)
        REGISTER(interaction, tag_detail::bosonic)
        
#undef REGISTER
        /**********************************************************************/
        
        generate_terms();
    }
    
    void update(BaseParameters const& p)
    {
        model << p;
        generate_terms();
        return;
    }
    
    Index<U1> const& phys_dim(size_t type) const
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
    typename U1::charge total_quantum_numbers(BaseParameters & parms) const
    {
        return static_cast<int>(parms["N_total"]);
    }
    
    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "n")
            return count;
        else if (name == "bdag")
            return create;
        else if (name == "b")
            return destroy;
        else if (name == "id")
            return ident;
        else if (name == "fill")
            return ident;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return 0;
    }
    
    table_ptr operators_table() const
    {
        return tag_handler;
    }
    
    measurements_type measurements () const
    {
        typedef std::vector<op_t> op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
        
        std::size_t ntypes = lat.maximum_vertex_type()+1;
        
        measurements_type meas;
        
        if (model["MEASURE[Density]"]) {
            std::string name = "Density";
            meas.push_back( new measurements::average<Matrix, U1>(name, lat,
                                                                  op_vec(ntypes,this->identity_matrix(0)),
                                                                  op_vec(ntypes,this->filling_matrix(0)),
                                                                  op_vec(ntypes,tag_handler->get_op(count))) );
        }
        
        if (model["MEASURE[Local density]"]) {
            std::string name = "Local density";
            meas.push_back( new measurements::local<Matrix, U1>(name, lat,
                                                                op_vec(ntypes,this->identity_matrix(0)),
                                                                op_vec(ntypes,this->filling_matrix(0)),
                                                                op_vec(ntypes,tag_handler->get_op(count))) );
        }
        
        if (model["MEASURE[Onebody density matrix]"]) {
            std::string name = "Onebody density matrix";
            bond_element ops;
            ops.push_back( std::make_pair(op_vec(ntypes,tag_handler->get_op(create)), false) );
            ops.push_back( std::make_pair(op_vec(ntypes,tag_handler->get_op(destroy)), false) );
            meas.push_back( new measurements::correlations<Matrix, U1>(name, lat,
                                                                       op_vec(ntypes,this->identity_matrix(0)),
                                                                       op_vec(ntypes,this->filling_matrix(0)),
                                                                       ops, true, false) );
        }
        
        return meas;
    }
    
private:
    void generate_terms() {
        this->terms_.clear();
        
        double k  = model["k"];
        double V0 = model["V0"];
        double w  = model["omega"];
        double shift = model["shift"];
        

        for (int p=0; p<lat.size(); ++p) {
            std::vector<int> neighs = lat.all(p);
            
            double x = lat.get_prop<double>("x", p);
            double exp_potential = V0*std::pow( std::cos(k*x), 2 );
            exp_potential += w*w/2. * std::pow(x - shift, 2 );
            
            double dx1 = lat.get_prop<double>("dx", p, neighs[0]);
            double dx2;
            if (neighs.size() == 1 && lat.get_prop<bool>("at_open_left_boundary", p))
                dx2 = lat.get_prop<double>("dx", p, p-1);
            else if (neighs.size() == 1 && lat.get_prop<bool>("at_open_right_boundary", p))
                dx2 = lat.get_prop<double>("dx", p, p+1);
            else
                dx2 = lat.get_prop<double>("dx", p, neighs[1]);
            
            double dx0 = lat.get_prop<double>("dx", p);
            
            // Psi''(x) = coeff1 * Psi(x+dx1) + coeff0 * Psi(x) + coeff2 * Psi(x+dx2)
            double coeff1 = 2. / (dx1*dx1 - dx1*dx2);
            double coeff2 = 2. / (dx2*dx2 - dx1*dx2);
            double coeff0 = -(coeff1 + coeff2);
            
            
            double U = model["c"] / dx0;
            double mu = exp_potential - model["mu"];
            mu += -coeff0 * model["h"];
            
            
#ifndef NDEBUG
            maquis::cout << "U = " << U << ", mu = " << mu << ", t = " << coeff1 * model["h"] << std::endl;
#endif
            
            if (U != 0.)
            { // interaction
                term_descriptor term;
                term.coeff = U;
                term.push_back( boost::make_tuple(p, interaction) );
                this->terms_.push_back(term);
            }
            
            if (mu != 0.)
            { // chemical potential
                term_descriptor term;
                term.coeff = mu;
                term.push_back( boost::make_tuple(p, count) );
                this->terms_.push_back(term);
            }
            
            for (int n=0; n<neighs.size(); ++n) { // hopping
                
                double t;
                if (lat.get_prop<double>("dx", p, neighs[n]) == dx1)
                    t = coeff1 * model["h"];
                else
                    t = coeff2 * model["h"];
                {
                    term_descriptor term;
                    term.coeff = -t;
                    term.push_back( boost::make_tuple(p, create) );
                    term.push_back( boost::make_tuple(neighs[n], destroy) );
                    this->terms_.push_back(term);
                }
            }
        }
    }

    
    const Lattice & lat;
    BaseParameters & model;
    Index<U1> phys;
    
    boost::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    tag_type ident, create, destroy, count, interaction;
};

#endif
