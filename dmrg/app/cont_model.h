//
//  model.h
//  
//
//  Created by Michele Dolfi on 30.06.11.
//  Copyright 2011 ETHZ. All rights reserved.
//

#ifndef CONTINUOUS_MODELS_H
#define CONTINUOUS_MODELS_H

#include "hamiltonian.h"
#include "measurements.h"

#include <cmath>

#include "utils/BaseParameters.h"

namespace app {
    
    /* ****************** OPTICAL LATTICE (with symmetry) */
    template<class Matrix>
    class OpticalLattice {
        typedef Hamiltonian<Matrix, U1> ham;        
        typedef typename ham::hamterm_t hamterm_t;        
        typedef Measurement_Term<Matrix, U1> mterm_t;
        typedef typename ham::op_t op_t;
    public:
        OpticalLattice (const Lattice& lat_, BaseParameters & model_)
        : lat(lat_)
        , model(model_)
        {
            phys.insert(std::make_pair(0, 1));
            ident.insert_block(Matrix(1, 1, 1), 0, 0);
            
            for (int n=1; n<=model.get<int>("Nmax"); ++n)
            {
                phys.insert(std::make_pair(n, 1));
                
                ident.insert_block(Matrix(1, 1, 1), n, n);
                
                count.insert_block(Matrix(1, 1, n), n, n);
                if ((n*n-n) != 0)
                    interaction.insert_block(Matrix(1, 1, n*n-n), n, n);
                
                
                create.insert_block(Matrix(1, 1, std::sqrt(n)), n-1, n);
                destroy.insert_block(Matrix(1, 1, std::sqrt(n)), n, n-1);
            }
        }
        
        Hamiltonian<Matrix, U1> H () const
        {
            std::vector<hamterm_t> terms;
            for (int p=0; p<lat.size(); ++p)
            {
                double exp_potential = model.get<double>("V0")*std::pow( std::sin(model.get<double>("k")*lat.get_prop<double>("x", p)), 2 );
                
                double U = model.get<double>("c")/lat.get_prop<double>("dx", p);
                double mu = ( -model.get<double>("mu")
                             + exp_potential
                             + model.get<double>("h")/(lat.get_prop<double>("dx", p)*lat.get_prop<double>("dx", p)) );
                if (!lat.get_prop<bool>("at_open_boundary", p))
                    mu += model.get<double>("h")/(lat.get_prop<double>("dx", p)*lat.get_prop<double>("dx", p));

                op_t site_op;
                for (int n=1; n<=model.get<int>("Nmax"); ++n)
                    if (U*(n*n-n)+mu*n != 0)
                        site_op.insert_block(Matrix(1, 1, U*(n*n-n)+mu*n), n, n);

                { // site term
                    hamterm_t term;
                    term.fill_operator = ident;
                    term.operators.push_back( std::make_pair(p, site_op) );
                    terms.push_back(term);
                }
                
//                if (mu != 0)
//                { // density term
//                    hamterm_t term;
//                    term.fill_operator = ident;
//                    term.operators.push_back( std::make_pair(p, mu*count) );
//                    terms.push_back(term);
//                }
//                
//                if (U != 0)
//                { // interaction term
//                    hamterm_t term;
//                    term.fill_operator = ident;
//                    term.operators.push_back( std::make_pair(p, U*interaction) );
//                    terms.push_back(term);
//                }
                
                std::vector<int> neighs = lat.forward(p);
                for (int n=0; n<neighs.size(); ++n) { // hopping
                    
                    double t = model.get<double>("h") / (lat.get_prop<double>("dx", p, neighs[n])*lat.get_prop<double>("dx", p, neighs[n]));
                    
                    {
                        hamterm_t term;
                        term.fill_operator = ident;
                        term.operators.push_back( std::make_pair(p, -t*create) );
                        term.operators.push_back( std::make_pair(neighs[n], destroy) );
                        terms.push_back(term);
                    }
                    {
                        hamterm_t term;
                        term.fill_operator = ident;
                        term.operators.push_back( std::make_pair(p, -t*destroy) );
                        term.operators.push_back( std::make_pair(neighs[n], create) );
                        terms.push_back(term);
                    }
                }
            }
            
            return ham(phys, ident, terms);
        }
        
        Measurements<Matrix, U1> measurements () const
        {
            Measurements<Matrix, U1> meas;
            meas.set_identity(ident);

            {
                mterm_t term;
                term.fill_operator = ident;
                term.name = "Density";
                term.type = mterm_t::Average;
                term.operators.push_back( std::make_pair(count, false) );
                
                meas.add_term(term);
            }
            {
                mterm_t term;
                term.fill_operator = ident;
                term.name = "Local density";
                term.type = mterm_t::Local;
                term.operators.push_back( std::make_pair(count, false) );
                
                meas.add_term(term);
            }
            return meas;
        }
        
    private:
        const Lattice & lat;
        BaseParameters & model;
        
        op_t ident;
        op_t create, destroy;
        op_t count, interaction;
        Index<U1> phys;
    };
    
    
    /* ****************** OPTICAL LATTICE (without symmetry) */
    template<class Matrix>
    class OpticalLatticeNull {
        typedef Hamiltonian<Matrix, NullGroup> ham;        
        typedef typename ham::hamterm_t hamterm_t;        
        typedef Measurement_Term<Matrix, NullGroup> mterm_t;
        typedef typename ham::op_t op_t;
    public:
        OpticalLatticeNull (const Lattice& lat_, BaseParameters & model_)
        : lat(lat_)
        , model(model_)
        {
            int phys_size = model.get<int>("Nmax")+1;
            NullGroup::charge c = NullGroup::SingletCharge;
            
            phys.insert(std::make_pair(c, phys_size));
            ident.insert_block(Matrix::identity_matrix(phys_size), c, c);
            count.insert_block(Matrix(phys_size, phys_size, 0), c, c);
            interaction.insert_block(Matrix(phys_size, phys_size, c), c, c);
            create.insert_block(Matrix(phys_size, phys_size, 0), c, c);
            destroy.insert_block(Matrix(phys_size, phys_size, 0), c, c);
            
            for (int n=1; n<=model.get<int>("Nmax"); ++n)
            {                                
                count[0](n, n) = n;

                interaction[0](n, n) = n*n-n;

                create[0](n-1, n) = std::sqrt(n);
                destroy[0](n, n-1) = std::sqrt(n);
            }
        }
        
        Hamiltonian<Matrix, NullGroup> H () const
        {
            NullGroup::charge c = NullGroup::SingletCharge;

            std::vector<hamterm_t> terms;
            for (int p=0; p<lat.size(); ++p)
            {
                double exp_potential = model.get<double>("V0")*std::pow( std::sin(model.get<double>("k")*lat.get_prop<double>("x", p)), 2 );
                
                double U = model.get<double>("c")/lat.get_prop<double>("dx", p);
                double mu = ( -model.get<double>("mu")
                             + exp_potential
                             + model.get<double>("h")/(lat.get_prop<double>("dx", p)*lat.get_prop<double>("dx", p)) );
                if (!lat.get_prop<bool>("at_open_boundary", p))
                    mu += model.get<double>("h")/(lat.get_prop<double>("dx", p)*lat.get_prop<double>("dx", p));
                
                op_t site_op;
                site_op.insert_block(U*interaction[0]+mu*count[0], c, c);
                
                { // site term
                    hamterm_t term;
                    term.fill_operator = ident;
                    term.operators.push_back( std::make_pair(p, site_op) );
                    terms.push_back(term);
                }
                
                //                if (mu != 0)
                //                { // density term
                //                    hamterm_t term;
                //                    term.fill_operator = ident;
                //                    term.operators.push_back( std::make_pair(p, mu*count) );
                //                    terms.push_back(term);
                //                }
                //                
                //                if (U != 0)
                //                { // interaction term
                //                    hamterm_t term;
                //                    term.fill_operator = ident;
                //                    term.operators.push_back( std::make_pair(p, U*interaction) );
                //                    terms.push_back(term);
                //                }
                
                std::vector<int> neighs = lat.forward(p);
                for (int n=0; n<neighs.size(); ++n) { // hopping
                    
                    double t = model.get<double>("h") / (lat.get_prop<double>("dx", p, neighs[n])*lat.get_prop<double>("dx", p, neighs[n]));
                    
                    {
                        hamterm_t term;
                        term.fill_operator = ident;
                        term.operators.push_back( std::make_pair(p, -t*create) );
                        term.operators.push_back( std::make_pair(neighs[n], destroy) );
                        terms.push_back(term);
                    }
                    {
                        hamterm_t term;
                        term.fill_operator = ident;
                        term.operators.push_back( std::make_pair(p, -t*destroy) );
                        term.operators.push_back( std::make_pair(neighs[n], create) );
                        terms.push_back(term);
                    }
                }
            }
            
            return ham(phys, ident, terms);
        }
        
        Measurements<Matrix, NullGroup> measurements () const
        {
            Measurements<Matrix, NullGroup> meas;
            meas.set_identity(ident);
            
            {
                mterm_t term;
                term.fill_operator = ident;
                term.name = "Density";
                term.type = mterm_t::Average;
                term.operators.push_back( std::make_pair(count, false) );
                
                meas.add_term(term);
            }
            {
                mterm_t term;
                term.fill_operator = ident;
                term.name = "Local density";
                term.type = mterm_t::Local;
                term.operators.push_back( std::make_pair(count, false) );
                
                meas.add_term(term);
            }
            return meas;
        }
        
    private:
        const Lattice & lat;
        BaseParameters & model;
        
        op_t ident;
        op_t create, destroy;
        op_t count, interaction;
        Index<NullGroup> phys;
    };

    
} // namespace

#endif
