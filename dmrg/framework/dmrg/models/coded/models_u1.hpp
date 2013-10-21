/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MODELS_CODED_U1_H
#define MODELS_CODED_U1_H

#include <sstream>

#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"

/* ****************** HEISENBERG */
template<class Matrix>
class Heisenberg : public Model<Matrix, U1>
{
    typedef Hamiltonian<Matrix, U1> ham;        
    typedef typename ham::hamterm_t hamterm_t;
    typedef typename ham::hamtagterm_t hamtagterm_t;
    typedef Measurement_Term<Matrix, U1> mterm_t;
    typedef typename ham::op_t op_t;
    typedef typename ham::table_type table_type;
    typedef typename ham::table_ptr table_ptr;
    typedef typename table_type::tag_type tag_type;
    typedef typename Matrix::value_type value_type;

public:   
    Heisenberg (const Lattice& lat, double Jxy, double Jz)
    : tag_handler(new table_type())
    {
        phys.insert(std::make_pair(1, 1));
        phys.insert(std::make_pair(-1, 1));

        op_t ident_op, splus_op, sminus_op, sz_op;
   
        ident_op.insert_block(Matrix(1, 1, 1), -1, -1);
        ident_op.insert_block(Matrix(1, 1, 1), 1, 1);
        
        splus_op.insert_block(Matrix(1, 1, 1), -1, 1);
        
        sminus_op.insert_block(Matrix(1, 1, 1), 1, -1);
        
        sz_op.insert_block(Matrix(1, 1, 0.5), 1, 1);
        sz_op.insert_block(Matrix(1, 1, -0.5), -1, -1);
        
        /**********************************************************************/
        /*** Create operator tag table ****************************************/
        /**********************************************************************/
        
#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);
        
        REGISTER(ident,   tag_detail::bosonic)
        REGISTER(splus,   tag_detail::bosonic)
        REGISTER(sminus,  tag_detail::bosonic)
        REGISTER(sz,      tag_detail::bosonic)
        
#undef REGISTER
        /**********************************************************************/

        
        for (int p=0; p<lat.size(); ++p) {
            std::vector<int> neighs = lat.forward(p);
            for (int n=0; n<neighs.size(); ++n) {
                {
                    hamtagterm_t term;
                    term.fill_operator = ident;
                    term.scale = Jz;
                    term.operators.push_back( std::make_pair(p, sz) );
                    term.operators.push_back( std::make_pair(neighs[n], sz) );
                    terms.push_back(term);
                }
                {
                    hamtagterm_t term;
                    term.fill_operator = ident;
                    term.scale = Jxy/2;
                    term.operators.push_back( std::make_pair(p, splus) );
                    term.operators.push_back( std::make_pair(neighs[n], sminus) );
                    terms.push_back(term);
                }
                {
                    hamtagterm_t term;
                    term.fill_operator = ident;
                    term.scale = Jxy/2;
                    term.operators.push_back( std::make_pair(p, sminus) );
                    term.operators.push_back( std::make_pair(neighs[n], splus) );
                    terms.push_back(term);
                }
            }
        }
        
    }
    
    Index<U1> get_phys() const
    {
        return phys;
    }

    Hamiltonian<Matrix, U1> H () const
    {        
        std::vector<hamterm_t> terms_ops;
        return ham(phys, tag_handler->get_op(ident), terms_ops, ident, terms, tag_handler);
    }
    
    Measurements<Matrix, U1> measurements () const
    {
        return Measurements<Matrix, U1>();
    }

    op_t get_op(std::string const & op) const
    {
        if (op == "splus")
            return tag_handler->get_op(splus);
        else if (op == "sminus")
            return tag_handler->get_op(sminus);
        else if (op == "sz")
            return tag_handler->get_op(sz);
        else
            throw std::runtime_error("Operator not valid for this model.");
        return op_t();
    }

private:
    Index<U1> phys;

    boost::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    tag_type ident, splus, sminus, sz;

    std::vector<hamtagterm_t> terms;
};

/* ****************** HARD CORE BOSONS */
template<class Matrix>
class HCB : public Model<Matrix, U1>
{
    typedef Hamiltonian<Matrix, U1> ham;
    typedef typename ham::hamterm_t hamterm_t;
    typedef typename ham::hamtagterm_t hamtagterm_t;
    typedef Measurement_Term<Matrix, U1> mterm_t;
    typedef typename ham::op_t op_t;
    typedef typename ham::table_type table_type;
    typedef typename ham::table_ptr table_ptr;
    typedef typename table_type::tag_type tag_type;
    typedef typename Matrix::value_type value_type;
    
public:   
    HCB (const Lattice& lat, double t=1)
    : tag_handler(new table_type())
    {
        phys.insert(std::make_pair(0, 1));
        phys.insert(std::make_pair(1, 1));

        op_t ident_op;
        op_t create_op, destroy_op, count_op;
        
        ident_op.insert_block(Matrix(1, 1, 1), 0, 0);
        ident_op.insert_block(Matrix(1, 1, 1), 1, 1);
        
        create_op.insert_block(Matrix(1, 1, 1), 0, 1);
        destroy_op.insert_block(Matrix(1, 1, 1), 1, 0);
        
        count_op.insert_block(Matrix(1, 1, 1), 1, 1);
        
        /**********************************************************************/
        /*** Create operator tag table ****************************************/
        /**********************************************************************/
        
#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);
        
        REGISTER(ident,   tag_detail::bosonic)
        REGISTER(create,  tag_detail::bosonic)
        REGISTER(destroy, tag_detail::bosonic)
        REGISTER(count,   tag_detail::bosonic)
        
#undef REGISTER
        /**********************************************************************/

        
        for (int p=0; p<lat.size(); ++p) {
            std::vector<int> neighs = lat.forward(p);
            for (int n=0; n<neighs.size(); ++n) {
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
            }
        }
        
    }
    
    
    Index<U1> get_phys() const
    {
        return phys;
    }

    Hamiltonian<Matrix, U1> H () const
    {        
        std::vector<hamterm_t> terms_ops;
        return ham(phys, tag_handler->get_op(ident), terms_ops, ident, terms, tag_handler);
    }
    
    Measurements<Matrix, U1> measurements () const
    {
        return Measurements<Matrix, U1>();
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
    Index<U1> phys;
    
    boost::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    tag_type ident, create, destroy, count;

    std::vector<hamtagterm_t> terms;
};

/* ****************** BOSE-HUBBARD */
template<class Matrix>
class BoseHubbard : public Model<Matrix, U1>
{
    typedef Hamiltonian<Matrix, U1> ham;
    typedef typename ham::hamterm_t hamterm_t;
    typedef typename ham::hamtagterm_t hamtagterm_t;
    typedef Measurement_Term<Matrix, U1> mterm_t;
    typedef typename ham::op_t op_t;
    typedef typename ham::table_type table_type;
    typedef typename ham::table_ptr table_ptr;
    typedef typename table_type::tag_type tag_type;
    typedef typename Matrix::value_type value_type;
    
public:
    BoseHubbard (const Lattice& lat_, BaseParameters & model_)
    : lat(lat_)
    , model(model_)
    , tag_handler(new table_type())
    {
        int Nmax = model["Nmax"];
        double t = model["t"];
        double U = model["U"];
        double V = model["V"];
        
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
            
            
            create_op.insert_block(Matrix(1, 1, std::sqrt(n)), n-1, n);
            destroy_op.insert_block(Matrix(1, 1, std::sqrt(n)), n, n-1);
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
    
    
    Index<U1> get_phys() const
    {
        return phys;
    }
    
    Hamiltonian<Matrix, U1> H () const
    {
        std::vector<hamterm_t> terms_ops;
        return ham(phys, tag_handler->get_op(ident), terms_ops, ident, terms, tag_handler);
    }
    
    Measurements<Matrix, U1> measurements () const
    {
        Measurements<Matrix, U1> meas;
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
        else if (op == "id")
            return tag_handler->get_op(ident);
        else if (op == "fill")
            return tag_handler->get_op(ident);
        else
            throw std::runtime_error("Operator not valid for this model.");
        return op_t();
    }
    
    
private:
    const Lattice & lat;
    BaseParameters & model;
    Index<U1> phys;
    
    boost::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    tag_type ident, create, destroy, count, interaction;
    
    std::vector<hamtagterm_t> terms;
};

/* ****************** FREE FERMIONS */
template<class Matrix>
class FreeFermions : public Model<Matrix, U1>
{
    typedef Hamiltonian<Matrix, U1> ham;
    typedef typename ham::hamterm_t hamterm_t;
    typedef typename ham::hamtagterm_t hamtagterm_t;
    typedef Measurement_Term<Matrix, U1> mterm_t;
    typedef typename ham::op_t op_t;
    typedef typename ham::table_type table_type;
    typedef typename ham::table_ptr table_ptr;
    typedef typename table_type::tag_type tag_type;
    typedef typename Matrix::value_type value_type;
    
public:
    FreeFermions (const Lattice& lat, double t=1)
    {
        op_t ident_op;
        op_t create_op, destroy_op, sign_op, dens_op;

        create_op.insert_block(Matrix(1, 1, 1), 0, 1);
        destroy_op.insert_block(Matrix(1, 1, 1), 1, 0);
        
        dens_op.insert_block(Matrix(1, 1, 1), 1, 1);
        
        ident_op.insert_block(Matrix(1, 1, 1), 0, 0);
        ident_op.insert_block(Matrix(1, 1, 1), 1, 1);
        sign_op.insert_block(Matrix(1, 1, 1), 0, 0);
        sign_op.insert_block(Matrix(1, 1, -1), 1, 1);
        
        phys.insert(std::make_pair(0, 1));
        phys.insert(std::make_pair(1, 1));
        
        /**********************************************************************/
        /*** Create operator tag table ****************************************/
        /**********************************************************************/
        
#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);
        
        REGISTER(ident,     tag_detail::bosonic)
        REGISTER(create,    tag_detail::fermionic)
        REGISTER(destroy,   tag_detail::fermionic)
        REGISTER(dens,      tag_detail::bosonic)
        REGISTER(sign,      tag_detail::bosonic)
        
#undef REGISTER
        /**********************************************************************/

        
        for (int p=0; p<lat.size(); ++p) {
            std::vector<int> neighs = lat.forward(p);
            for (int n=0; n<neighs.size(); ++n) {
                {
                    hamtagterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign;
                    term.scale = -t;
                    term.operators.push_back( std::make_pair(p, create) );
                    term.operators.push_back( std::make_pair(neighs[n], destroy) );
                    terms.push_back(term);
                }
                {
                    hamtagterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign;
                    term.scale = -t;
                    term.operators.push_back( std::make_pair(p, destroy) );
                    term.operators.push_back( std::make_pair(neighs[n], create) );
                    terms.push_back(term);
                }
            }
        }
        
    }
    
    Index<U1> get_phys() const
    {
        return phys;
    }

    Hamiltonian<Matrix, U1> H () const
    {        
        std::vector<hamterm_t> terms_ops;
        return ham(phys, tag_handler->get_op(ident), terms_ops, ident, terms, tag_handler);
    }
    
    Measurements<Matrix, U1> measurements () const
    {
        Measurements<Matrix, U1> meas;
        meas.set_identity(tag_handler->get_op(ident));
        
        {
            mterm_t term;
            term.name = "Density";
            term.type = mterm_t::Local;
            term.operators.push_back( std::make_pair(tag_handler->get_op(dens), false) );
            
            meas.add_term(term);
        }
        {
            mterm_t term;
            term.name = "DensityCorrelation";
            term.type = mterm_t::HalfCorrelation;
            term.fill_operator = tag_handler->get_op(ident);
            term.operators.push_back( std::make_pair(tag_handler->get_op(dens), false) );
            term.operators.push_back( std::make_pair(tag_handler->get_op(dens), false) );
            
            meas.add_term(term);
        }
        {
            mterm_t term;
            term.name = "OneBodyDM";
            term.type = mterm_t::HalfCorrelation;
            term.fill_operator = tag_handler->get_op(sign);
            term.operators.push_back( std::make_pair(tag_handler->get_op(create), true) );
            term.operators.push_back( std::make_pair(tag_handler->get_op(destroy), true) );
            
            meas.add_term(term);
        }
        
        return meas;
    }
    
    op_t get_op(std::string const & op) const
    {
        if (op == "n")
            return tag_handler->get_op(dens);
        else if (op == "cdag")
            return tag_handler->get_op(create);
        else if (op == "c")
            return tag_handler->get_op(destroy);
        else
            throw std::runtime_error("Operator not valid for this model.");
        return op_t();
    }
    
private:
    Index<U1> phys;

    boost::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    tag_type ident;
    tag_type create, destroy, sign, dens;

    std::vector<hamtagterm_t> terms;
};


/* ****************** FERMI HUBBARD */
template<class Matrix>
class FermiHubbardU1 : public Model<Matrix, U1>
{
public:
    typedef Hamiltonian<Matrix, U1> ham;
    typedef typename ham::hamterm_t hamterm_t;
    typedef typename ham::op_t op_t;
    typedef typename ham::hamtagterm_t hamtagterm_t;
    typedef Measurement_Term<Matrix, U1> mterm_t;
    typedef typename ham::table_type table_type;
    typedef typename ham::table_ptr table_ptr;
    typedef typename table_type::tag_type tag_type;
    typedef typename Matrix::value_type value_type;
    
    FermiHubbardU1(const Lattice& lat, BaseParameters & parms)
    {
        phys.insert(std::make_pair(0, 1));
        phys.insert(std::make_pair(1, 2));
        phys.insert(std::make_pair(2, 1));
        op_t create_up_op, create_down_op, destroy_up_op, destroy_down_op,
            count_up_op, count_down_op, doubly_occ_op,
            sign_up_op, sign_down_op, fill_op, ident_op;
        
        ident_op.insert_block(Matrix(1, 1, 1), 0, 0);
        ident_op.insert_block(Matrix::identity_matrix(2), 1, 1);
        ident_op.insert_block(Matrix(1, 1, 1), 2, 2);
        
        {Matrix tmp(1,2,0); tmp(0,0)=1; create_up_op.insert_block(tmp, 0, 1);}
        {Matrix tmp(2,1,0); tmp(1,0)=1; create_up_op.insert_block(tmp, 1, 2);}
        {Matrix tmp(1,2,0); tmp(0,1)=1; create_down_op.insert_block(tmp, 0, 1);}
        {Matrix tmp(2,1,0); tmp(0,0)=1; create_down_op.insert_block(tmp, 1, 2);}
        
        {Matrix tmp(2,1,0); tmp(0,0)=1; destroy_up_op.insert_block(tmp, 1, 0);}
        {Matrix tmp(1,2,0); tmp(0,1)=1; destroy_up_op.insert_block(tmp, 2, 1);}
        {Matrix tmp(2,1,0); tmp(1,0)=1; destroy_down_op.insert_block(tmp, 1, 0);}
        {Matrix tmp(1,2,0); tmp(0,0)=1; destroy_down_op.insert_block(tmp, 2, 1);}
        
        {Matrix tmp(2,2,0); tmp(0,0)=1; count_up_op.insert_block(tmp, 1, 1);}
        count_up_op.insert_block(Matrix(1, 1, 1), 2, 2);
        {Matrix tmp(2,2,0); tmp(1,1)=1; count_down_op.insert_block(tmp, 1, 1);}
        count_down_op.insert_block(Matrix(1, 1, 1), 2, 2);
        
        doubly_occ_op.insert_block(Matrix(1, 1, 1), 2, 2);
        
        sign_up_op.insert_block(Matrix(1, 1, 1), 0, 0);
        {Matrix tmp=Matrix::identity_matrix(2); tmp(0,0)=-1; sign_up_op.insert_block(tmp, 1, 1);}
        sign_up_op.insert_block(Matrix(1, 1, -1), 2, 2);
        
        sign_down_op.insert_block(Matrix(1, 1, 1), 0, 0);
        {Matrix tmp=Matrix::identity_matrix(2); tmp(1,1)=-1; sign_down_op.insert_block(tmp, 1, 1);}
        sign_down_op.insert_block(Matrix(1, 1, -1), 2, 2);
        
        gemm(sign_up_op, sign_down_op, fill_op);

        /**********************************************************************/
        /*** Create operator tag table ****************************************/
        /**********************************************************************/
        
#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);
        
        REGISTER(ident,         tag_detail::bosonic)
        REGISTER(fill,          tag_detail::bosonic)
        REGISTER(create_up,     tag_detail::fermionic)
        REGISTER(create_down,   tag_detail::fermionic)
        REGISTER(destroy_up,    tag_detail::fermionic)
        REGISTER(destroy_down,  tag_detail::fermionic)
        REGISTER(count_up,      tag_detail::bosonic)
        REGISTER(count_down,    tag_detail::bosonic)
        REGISTER(doubly_occ,    tag_detail::bosonic)
        REGISTER(sign_up,       tag_detail::bosonic)
        REGISTER(sign_down,     tag_detail::bosonic)
        
#undef REGISTER
        /**********************************************************************/
        
        double U = parms["U"];
        op_t tmp;
        for (int p=0; p<lat.size(); ++p) {
            { // U term
                hamtagterm_t term;
                term.fill_operator = ident;
                term.scale = U;
                term.operators.push_back( std::make_pair(p, doubly_occ) );
                terms.push_back(term);
            }
            
            std::vector<int> neighs = lat.forward(p);
            for (std::vector<int>::iterator hopto = neighs.begin();
                 hopto != neighs.end(); ++hopto)
            {
                double ti = get_t(parms,
                                  lat.get_prop<int>("type", p, *hopto));
                { // t*cdag_up*c_up
                    hamtagterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_up;

                    // Note inverse notation because of notation in operator.
                    std::pair<tag_type, value_type> prod = tag_handler->get_product_tag(sign_up, create_up);
                    term.scale = -ti * prod.second;
                    term.operators.push_back( std::make_pair(p, prod.first) );
                    term.operators.push_back( std::make_pair(*hopto, destroy_up) );
                    terms.push_back(term);
                }
                { // t*c_up*cdag_up
                    hamtagterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_up;

                    // Note inverse notation because of notation in operator.
                    std::pair<tag_type, value_type> prod = tag_handler->get_product_tag(destroy_up, sign_up);
                    term.scale = -ti * prod.second;
                    term.operators.push_back( std::make_pair(p, prod.first) );
                    term.operators.push_back( std::make_pair(*hopto, create_up) );
                    terms.push_back(term);
                }
                { // t*cdag_down*c_down
                    hamtagterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_down;

                    // Note inverse notation because of notation in operator.
                    std::pair<tag_type, value_type> prod = tag_handler->get_product_tag(sign_down, create_down);
                    term.scale = -ti * prod.second;
                    term.operators.push_back( std::make_pair(p, prod.first) );
                    term.operators.push_back( std::make_pair(*hopto, destroy_down) );
                    terms.push_back(term);
                }
                { // t*c_down*cdag_down
                    hamtagterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_down;

                    // Note inverse notation because of notation in operator.
                    std::pair<tag_type, value_type> prod = tag_handler->get_product_tag(destroy_down, sign_down);
                    term.scale = -ti * prod.second;
                    term.operators.push_back( std::make_pair(p, prod.first) );
                    term.operators.push_back( std::make_pair(*hopto, create_down) );
                    terms.push_back(term);
                }
            }
        }
    }
    
    Index<U1> get_phys() const
    {
        return phys;
    }
    
    Hamiltonian<Matrix, U1> H () const
    {
        std::vector<hamterm_t> terms_ops;
        return ham(phys, tag_handler->get_op(ident), terms_ops, ident, terms, tag_handler);
    }
    
    Measurements<Matrix, U1> measurements () const
    {
        return Measurements<Matrix, U1>();
    }

    op_t get_op(std::string const & op) const
    {
        if (op == "n_up")
            return tag_handler->get_op(count_up);
        else if (op == "n_down")
            return tag_handler->get_op(count_down);
        else if (op == "cdag_up")
            return tag_handler->get_op(create_up);
        else if (op == "cdag_down")
            return tag_handler->get_op(create_down);
        else if (op == "c_up")
            return tag_handler->get_op(destroy_up);
        else if (op == "c_down")
            return tag_handler->get_op(destroy_down);
        else if (op == "id")
            return tag_handler->get_op(ident);
        else if (op == "fill")
            return tag_handler->get_op(fill);
        else
            throw std::runtime_error("Operator not valid for this model.");
        return op_t();
    }

    
private:
    Index<U1> phys;
    tag_type create_up, create_down, destroy_up, destroy_down, count_up, count_down, doubly_occ,
             sign_up, sign_down, fill, ident;

    std::vector<hamtagterm_t> terms;

    boost::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    
    double get_t (BaseParameters & parms, int i)
    {
        std::ostringstream key;
        key << "t" << (i+1);
        return (parms.is_set(key.str())) ? parms[key.str()] : parms["t"];
    }
};

#endif
