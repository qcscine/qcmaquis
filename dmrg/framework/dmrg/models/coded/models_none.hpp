/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MODELS_CODED_NONE_H
#define MODELS_CODED_NONE_H

#include <sstream>

#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"

/* ****************** BOSE-HUBBARD */
template<class Matrix>
class BoseHubbardNone : public Model<Matrix, TrivialGroup>
{
    typedef Hamiltonian<Matrix, TrivialGroup> ham;
    typedef typename ham::hamterm_t hamterm_t;
    typedef typename ham::hamtagterm_t hamtagterm_t;
    typedef Measurement_Term<Matrix, TrivialGroup> mterm_t;
    typedef typename ham::op_t op_t;
    typedef typename ham::table_type table_type;
    typedef typename ham::table_ptr table_ptr;
    typedef typename table_type::tag_type tag_type;
    typedef typename Matrix::value_type value_type;
    
public:
    BoseHubbardNone (const Lattice& lat, BaseParameters & model_)
    : model(model_)
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
            
            mcreate(n-1,n) = std::sqrt(n);   // input n-1, output n
            mdestroy(n,n-1) = std::sqrt(n);  // input n,   output n-1
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
            {
                hamtagterm_t term;
                term.with_sign = false;
                term.fill_operator = ident;
                term.scale = U/2.;
                term.operators.push_back( std::make_pair(p, interaction) );
                terms.push_back(term);
            }
            
            std::vector<int> neighs = lat.forward(p);
            for (int n=0; n<neighs.size(); ++n) {
                /* hopping */
                {
                    hamtagterm_t term;
                    term.fill_operator = ident;
                    term.scale = -t;
                    term.operators.push_back( std::make_pair(p, create) );
                    term.operators.push_back( std::make_pair(neighs[n], destroy) );
                    terms.push_back(term);
                }
                {
                    hamtagterm_t term;
                    term.fill_operator = ident;
                    term.scale = -t;
                    term.operators.push_back( std::make_pair(p, destroy) );
                    term.operators.push_back( std::make_pair(neighs[n], create) );
                    terms.push_back(term);
                }
                /* nearest-neighborn interaction */
                {
                    hamtagterm_t term;
                    term.fill_operator = ident;
                    term.scale = V;
                    term.operators.push_back( std::make_pair(p, count) );
                    term.operators.push_back( std::make_pair(neighs[n], count) );
                    terms.push_back(term);
                }
            }
        }
        
    }
    
    Index<TrivialGroup> get_phys() const
    {
        return phys;
    }
    
    Hamiltonian<Matrix, TrivialGroup> H () const
    {
        std::vector<hamterm_t> terms_ops;
        return ham(phys, tag_handler->get_op(ident), terms_ops, ident, terms, tag_handler);
    }
    
    Measurements<Matrix, TrivialGroup> measurements () const
    {
        Measurements<Matrix, TrivialGroup> meas;
        meas.set_identity(tag_handler->get_op(ident));
        
        if (model["ENABLE_MEASURE[Density]"]) {
            mterm_t term;
            term.fill_operator = tag_handler->get_op(ident);
            term.name = "Density";
            term.type = mterm_t::Average;
            term.operators.push_back( std::make_pair(tag_handler->get_op(count), false) );
            
            meas.add_term(term);
        }
        if (model["ENABLE_MEASURE[Local density]"]) {
            mterm_t term;
            term.fill_operator = tag_handler->get_op(ident);
            term.name = "Local density";
            term.type = mterm_t::Local;
            term.operators.push_back( std::make_pair(tag_handler->get_op(count), false) );
            
            meas.add_term(term);
        }
        
        if (model["ENABLE_MEASURE[Onebody density matrix]"]) {
            mterm_t term;
            term.fill_operator = tag_handler->get_op(ident);
            term.name = "Onebody density matrix";
            term.type = mterm_t::HalfCorrelation;
            term.operators.push_back( std::make_pair(tag_handler->get_op(create), false) );
            term.operators.push_back( std::make_pair(tag_handler->get_op(destroy), false) );
            
            meas.add_term(term);
        }
        
        return meas;
    }
    
    op_t get_op(std::string const & op) const
    {
        if (op == "n")
            return tag_handler->get_op(count);
        else if (op == "bdag")
            return tag_handler->get_op(create);
        else if (op == "b")
            return tag_handler->get_op(destroy);
        else
            throw std::runtime_error("Operator not valid for this model.");
        return op_t();
    }
    
    
private:
    BaseParameters & model;
    Index<TrivialGroup> phys;

    boost::shared_ptr<TagHandler<Matrix, TrivialGroup> > tag_handler;
    tag_type ident, create, destroy, count, interaction;

    
    std::vector<hamtagterm_t> terms;
};



#endif
