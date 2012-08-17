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

#endif
