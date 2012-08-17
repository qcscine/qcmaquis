/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_HAMILTONIANS_H
#define APP_HAMILTONIANS_H

#include <sstream>

#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"

/* ****************** HEISENBERG */
template<class Matrix>
class Heisenberg : public Model<Matrix, U1>
{
    typedef Hamiltonian<Matrix, U1> ham;        
    typedef typename ham::hamterm_t hamterm_t;        
    typedef typename ham::op_t op_t;

public:   
    Heisenberg (const Lattice& lat, double Jxy, double Jz)
    {
        
        ident.insert_block(Matrix(1, 1, 1), -1, -1);
        ident.insert_block(Matrix(1, 1, 1), 1, 1);
        
        splus.insert_block(Matrix(1, 1, 1), -1, 1);
        
        sminus.insert_block(Matrix(1, 1, 1), 1, -1);
        
        sz.insert_block(Matrix(1, 1, 0.5), 1, 1);
        sz.insert_block(Matrix(1, 1, -0.5), -1, -1);
        
        phys.insert(std::make_pair(1, 1));
        phys.insert(std::make_pair(-1, 1));
        
        for (int p=0; p<lat.size(); ++p) {
            std::vector<int> neighs = lat.forward(p);
            for (int n=0; n<neighs.size(); ++n) {
                {
                    hamterm_t term;
                    term.fill_operator = ident;
                    term.operators.push_back( std::make_pair(p, Jz*sz) );
                    term.operators.push_back( std::make_pair(neighs[n], sz) );
                    terms.push_back(term);
                }
                {
                    hamterm_t term;
                    term.fill_operator = ident;
                    term.operators.push_back( std::make_pair(p, Jxy/2*splus) );
                    term.operators.push_back( std::make_pair(neighs[n], sminus) );
                    terms.push_back(term);
                }
                {
                    hamterm_t term;
                    term.fill_operator = ident;
                    term.operators.push_back( std::make_pair(p, Jxy/2*sminus) );
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
        return ham(phys, ident, terms);
    }
    
    Measurements<Matrix, U1> measurements () const
    {
        return Measurements<Matrix, U1>();
    }

    op_t get_op(std::string const & op) const
    {
        if (op == "splus")
            return splus;
        else if (op == "sminus")
            return sminus;
        else if (op == "sz")
            return sz;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return op_t();
    }

private:
    op_t ident;    
    op_t splus, sminus, sz;
    Index<U1> phys;

    std::vector<hamterm_t> terms;
};

/* ****************** HARD CORE BOSONS */
template<class Matrix>
class HCB : public Model<Matrix, U1>
{
    typedef Hamiltonian<Matrix, U1> ham;        
    typedef typename ham::hamterm_t hamterm_t;        
    typedef typename ham::op_t op_t;
    
public:   
    HCB (const Lattice& lat, double t=1)
    {
        ident.insert_block(Matrix(1, 1, 1), 0, 0);
        ident.insert_block(Matrix(1, 1, 1), 1, 1);
        
        create.insert_block(Matrix(1, 1, 1), 0, 1);
        destroy.insert_block(Matrix(1, 1, 1), 1, 0);
        
        count.insert_block(Matrix(1, 1, 1), 1, 1);
        
        phys.insert(std::make_pair(0, 1));
        phys.insert(std::make_pair(1, 1));
        
        for (int p=0; p<lat.size(); ++p) {
            std::vector<int> neighs = lat.forward(p);
            for (int n=0; n<neighs.size(); ++n) {
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
        
    }
    
    
    Index<U1> get_phys() const
    {
        return phys;
    }

    Hamiltonian<Matrix, U1> H () const
    {        
        return ham(phys, ident, terms);
    }
    
    Measurements<Matrix, U1> measurements () const
    {
        return Measurements<Matrix, U1>();
    }
    
    op_t get_op(std::string const & op) const
    {
        if (op == "n")
            return count;
        else if (op == "bdag")
            return create;
        else if (op == "b")
            return destroy;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return op_t();
    }
    
    
private:
    op_t ident;    
    op_t create, destroy, count;
    Index<U1> phys;
    
    std::vector<hamterm_t> terms;
};

/* ****************** BOSE-HUBBARD */
template<class Matrix>
class BoseHubbard : public Model<Matrix, U1>
{
    typedef Hamiltonian<Matrix, U1> ham;
    typedef typename ham::hamterm_t hamterm_t;
    typedef typename ham::op_t op_t;
    
public:
    BoseHubbard (const Lattice& lat, int Nmax=2, double t=1., double U=1., double V=1.)
    {
        phys.insert(std::make_pair(0, 1));
        ident.insert_block(Matrix(1, 1, 1), 0, 0);
        
        for (int n=1; n<=Nmax; ++n)
        {
            phys.insert(std::make_pair(n, 1));
            
            ident.insert_block(Matrix(1, 1, 1), n, n);
            
            count.insert_block(Matrix(1, 1, n), n, n);
            if ((n*n-n) != 0)
                interaction.insert_block(Matrix(1, 1, n*n-n), n, n);
            
            
            create.insert_block(Matrix(1, 1, std::sqrt(n)), n-1, n);
            destroy.insert_block(Matrix(1, 1, std::sqrt(n)), n, n-1);
        }
        
        for (int p=0; p<lat.size(); ++p) {
            /* interaction */
            {
                hamterm_t term;
                term.with_sign = false;
                term.fill_operator = ident;
                term.operators.push_back( std::make_pair(p, interaction) );
                terms.push_back(term);
            }
            
            std::vector<int> neighs = lat.forward(p);
            for (int n=0; n<neighs.size(); ++n) {
                /* hopping */
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
                /* nearest-neighborn interaction */
                {
                    hamterm_t term;
                    term.fill_operator = ident;
                    term.operators.push_back( std::make_pair(p, V*count) );
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
        return ham(phys, ident, terms);
    }
    
    Measurements<Matrix, U1> measurements () const
    {
        return Measurements<Matrix, U1>();
    }
    
    op_t get_op(std::string const & op) const
    {
        if (op == "n")
            return count;
        else if (op == "bdag")
            return create;
        else if (op == "b")
            return destroy;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return op_t();
    }
    
    
private:
    op_t ident;
    op_t create, destroy, count, interaction;
    Index<U1> phys;
    
    std::vector<hamterm_t> terms;
};

/* ****************** FREE FERMIONS */
template<class Matrix>
class FreeFermions : public Model<Matrix, U1>
{
    typedef Hamiltonian<Matrix, U1> ham;        
    typedef typename ham::hamterm_t hamterm_t;        
    typedef typename ham::op_t op_t;
    typedef Measurement_Term<Matrix, U1> mterm_t;
    
public:   
    FreeFermions (const Lattice& lat, double t=1)
{
    create.insert_block(Matrix(1, 1, 1), 0, 1);
    destroy.insert_block(Matrix(1, 1, 1), 1, 0);
    
    dens.insert_block(Matrix(1, 1, 1), 1, 1);
    
    ident.insert_block(Matrix(1, 1, 1), 0, 0);
    ident.insert_block(Matrix(1, 1, 1), 1, 1);
    sign.insert_block(Matrix(1, 1, 1), 0, 0);
    sign.insert_block(Matrix(1, 1, -1), 1, 1);
    
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    
    for (int p=0; p<lat.size(); ++p) {
        std::vector<int> neighs = lat.forward(p);
        for (int n=0; n<neighs.size(); ++n) {
            {
                hamterm_t term;
                term.with_sign = true;
                term.fill_operator = sign;
                term.operators.push_back( std::make_pair(p, -t*create) );
                term.operators.push_back( std::make_pair(neighs[n], destroy) );
                terms.push_back(term);
            }
            {
                hamterm_t term;
                term.with_sign = true;
                term.fill_operator = sign;
                term.operators.push_back( std::make_pair(p, -t*destroy) );
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
        return ham(phys, ident, terms);
    }
    
    Measurements<Matrix, U1> measurements () const
    {
        Measurements<Matrix, U1> meas;
        meas.set_identity(ident);
        
        {
            mterm_t term;
            term.name = "Density";
            term.type = mterm_t::Local;
            term.operators.push_back( std::make_pair(dens, false) );
            
            meas.add_term(term);
        }
        {
            mterm_t term;
            term.name = "DensityCorrelation";
            term.type = mterm_t::HalfCorrelation;
            term.fill_operator = ident;
            term.operators.push_back( std::make_pair(dens, false) );
            term.operators.push_back( std::make_pair(dens, false) );
            
            meas.add_term(term);
        }
        {
            mterm_t term;
            term.name = "OneBodyDM";
            term.type = mterm_t::HalfCorrelation;
            term.fill_operator = sign;
            term.operators.push_back( std::make_pair(create, true) );
            term.operators.push_back( std::make_pair(destroy, true) );
            
            meas.add_term(term);
        }
        
        return meas;
    }
    
    op_t get_op(std::string const & op) const
    {
        if (op == "n")
            return dens;
        else if (op == "cdag")
            return create;
        else if (op == "c")
            return destroy;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return op_t();
    }
    
private:
    op_t ident;    
    op_t create, destroy, sign, dens;
    Index<U1> phys;
    
    std::vector<hamterm_t> terms;
};


/* ****************** FERMI HUBBARD */
template<class Matrix>
class FermiHubbardU1 : public Model<Matrix, U1>
{
public:
    typedef Hamiltonian<Matrix, U1> ham;
    typedef typename ham::hamterm_t hamterm_t;
    typedef typename ham::op_t op_t;
    
    FermiHubbardU1(const Lattice& lat, BaseParameters & parms)
    {
        phys.insert(std::make_pair(0, 1));
        phys.insert(std::make_pair(1, 2));
        phys.insert(std::make_pair(2, 1));
        
        ident.insert_block(Matrix(1, 1, 1), 0, 0);
        ident.insert_block(identity_matrix(Matrix(),2), 1, 1);
        ident.insert_block(Matrix(1, 1, 1), 2, 2);
        
        {Matrix tmp(1,2,0); tmp(0,0)=1; create_up.insert_block(tmp, 0, 1);}
        {Matrix tmp(2,1,0); tmp(1,0)=1; create_up.insert_block(tmp, 1, 2);}
        {Matrix tmp(1,2,0); tmp(0,1)=1; create_down.insert_block(tmp, 0, 1);}
        {Matrix tmp(2,1,0); tmp(0,0)=1; create_down.insert_block(tmp, 1, 2);}
        
        {Matrix tmp(2,1,0); tmp(0,0)=1; destroy_up.insert_block(tmp, 1, 0);}
        {Matrix tmp(1,2,0); tmp(0,1)=1; destroy_up.insert_block(tmp, 2, 1);}
        {Matrix tmp(2,1,0); tmp(1,0)=1; destroy_down.insert_block(tmp, 1, 0);}
        {Matrix tmp(1,2,0); tmp(0,0)=1; destroy_down.insert_block(tmp, 2, 1);}
        
        {Matrix tmp(2,2,0); tmp(0,0)=1; count_up.insert_block(tmp, 1, 1);}
        count_up.insert_block(Matrix(1, 1, 1), 2, 2);
        {Matrix tmp(2,2,0); tmp(1,1)=1; count_down.insert_block(tmp, 1, 1);}
        count_down.insert_block(Matrix(1, 1, 1), 2, 2);
        
        doubly_occ.insert_block(Matrix(1, 1, 1), 2, 2);
        
        sign_up.insert_block(Matrix(1, 1, 1), 0, 0);
        {Matrix tmp=identity_matrix(Matrix(),2); tmp(0,0)=-1; sign_up.insert_block(tmp, 1, 1);}
        sign_up.insert_block(Matrix(1, 1, -1), 2, 2);
        
        sign_down.insert_block(Matrix(1, 1, 1), 0, 0);
        {Matrix tmp=identity_matrix(Matrix(),2); tmp(1,1)=-1; sign_down.insert_block(tmp, 1, 1);}
        sign_down.insert_block(Matrix(1, 1, -1), 2, 2);
        
        op_t tmp;
        for (int p=0; p<lat.size(); ++p) {
            { // U term
                hamterm_t term;
                term.fill_operator = ident;
                term.operators.push_back( std::make_pair(p, parms.get<double>("U")*doubly_occ) );
                terms.push_back(term);
            }
            
            std::vector<int> neighs = lat.forward(p);
            for (std::vector<int>::iterator hopto = neighs.begin();
                 hopto != neighs.end(); ++hopto)
            {
                double ti = get_t(parms,
                                  lat.get_prop<int>("type", p, *hopto));
                { // t*cdag_up*c_up
                    hamterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_up;
                    gemm(sign_up, create_up, tmp); // Note inverse notation because of notation in operator.
                    term.operators.push_back( std::make_pair(p, -ti*tmp) );
                    term.operators.push_back( std::make_pair(*hopto, destroy_up) );
                    terms.push_back(term);
                }
                { // t*c_up*cdag_up
                    hamterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_up;
                    gemm(destroy_up, sign_up, tmp); // Note inverse notation because of notation in operator.
                    term.operators.push_back( std::make_pair(p, -ti*tmp) );
                    term.operators.push_back( std::make_pair(*hopto, create_up) );
                    terms.push_back(term);
                }
                { // t*cdag_down*c_down
                    hamterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_down;
                    gemm(sign_down, create_down, tmp); // Note inverse notation because of notation in operator.
                    term.operators.push_back( std::make_pair(p, -ti*tmp) );
                    term.operators.push_back( std::make_pair(*hopto, destroy_down) );
                    terms.push_back(term);
                }
                { // t*c_down*cdag_down
                    hamterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_down;
                    gemm(destroy_down, sign_down, tmp); // Note inverse notation because of notation in operator.
                    term.operators.push_back( std::make_pair(p, -ti*tmp) );
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
        return ham(phys, ident, terms);
    }
    
    Measurements<Matrix, U1> measurements () const
    {
        return Measurements<Matrix, U1>();
    }

    op_t get_op(std::string const & op) const
    {
        if (op == "n_up")
            return count_up;
        else if (op == "n_down")
            return count_down;
        else if (op == "cdag_up")
            return create_up;
        else if (op == "cdag_down")
            return create_down;
        else if (op == "c_up")
            return destroy_up;
        else if (op == "c_down")
            return destroy_down;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return op_t();
    }

    
private:
    Index<U1> phys;
    op_t create_up, create_down, destroy_up, destroy_down;
    op_t count_up, count_down, doubly_occ;
    op_t sign_up, sign_down;
    op_t ident;
    std::vector<hamterm_t> terms;
    
    double get_t (BaseParameters & parms, int i)
    {
        std::ostringstream key;
        key << "t" << (i+1);
        return (parms.is_set(key.str())) ? parms.get<double>(key.str()) : parms.get<double>("t");
    }
};

#endif
