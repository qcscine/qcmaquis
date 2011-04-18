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

#include "hamiltonian.h"

namespace app {
    /* ****************** HEISENBERG */
    template<class Matrix>
    class Heisenberg : public Hamiltonian<Matrix, U1>
    {
    public:
        typedef Hamiltonian_Term<Matrix, U1> hamterm_t;
        typedef typename hamterm_t::op_t op_t;
        
        Heisenberg(const Lattice& lat, double Jxy_, double Jz_) : Jxy(Jxy_), Jz(Jz_)
        {
            op_t splus, sminus, sz;
            
            ident.insert_block(Matrix(1, 1, 1), -1, -1);
            ident.insert_block(Matrix(1, 1, 1), 1, 1);
            
            splus.insert_block(Matrix(1, 1, 1), -1, 1);
            
            sminus.insert_block(Matrix(1, 1, 1), 1, -1);
            
            sz.insert_block(Matrix(1, 1, 0.5), 1, 1);
            sz.insert_block(Matrix(1, 1, -0.5), -1, -1);
            
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
        
        int n_terms() const
        {
            return terms.size();
        }
        hamterm_t operator[](int i) const
        {
            return terms[i];
        }
        
        op_t get_identity() const
        {
            return ident;
        }
        
        Index<U1> get_phys() const
        {
            Index<U1> phys;
            phys.insert(std::make_pair(1, 1));
            phys.insert(std::make_pair(-1, 1));
            return phys;
        }
        
    private:
        double Jxy, Jz;
        std::vector<hamterm_t> terms;
        op_t ident;
    };
    
    
    /* ****************** HARD CORE BOSONS */
    template<class Matrix>
    class HCB : public Hamiltonian<Matrix, U1>
    {
    public:
        typedef Hamiltonian_Term<Matrix, U1> hamterm_t;
        typedef typename hamterm_t::op_t op_t;
        
        HCB (const Lattice& lat, double t=1.)
        {
            op_t create, destroy, count;
            
            ident.insert_block(Matrix(1, 1, 1), 0, 0);
            ident.insert_block(Matrix(1, 1, 1), 1, 1);
            
            create.insert_block(Matrix(1, 1, 1), 0, 1);
            destroy.insert_block(Matrix(1, 1, 1), 1, 0);
            
            count.insert_block(Matrix(1, 1, 1), 1, 1);
            
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
        
        int n_terms() const
        {
            return terms.size();
        }
        hamterm_t operator[](int i) const
        {
            return terms[i];
        }
        
        op_t get_identity() const
        {
            return ident;
        }
        
        Index<U1> get_phys() const
        {
            Index<U1> phys;
            phys.insert(std::make_pair(0, 1));
            phys.insert(std::make_pair(1, 1));
            return phys;
        }
        
    private:
        std::vector<hamterm_t> terms;
        op_t ident;
    };
    
    
    /* ****************** FERMI HUBBARD */
    template<class Matrix>
    class TwoU1_FermiHubbard : public Hamiltonian<Matrix, TwoU1>
    {
    public:
        typedef Hamiltonian_Term<Matrix, TwoU1> hamterm_t;
        typedef typename hamterm_t::op_t op_t;
        
        TwoU1_FermiHubbard(const Lattice& lat, BaseParameters & parms)
        {
            op_t create_up, create_down, destroy_up, destroy_down;
            op_t count_up, count_down, doubly_occ;
            op_t sign_up, sign_down;
            
            TwoU1::charge A(0), B(0), C(0), D(1);
            B[0]=1; C[1]=1;
            
            ident.insert_block(Matrix(1, 1, 1), A, A);
            ident.insert_block(Matrix(1, 1, 1), B, B);
            ident.insert_block(Matrix(1, 1, 1), C, C);
            ident.insert_block(Matrix(1, 1, 1), D, D);
            
            create_up.insert_block(Matrix(1, 1, 1), A, B);
            create_up.insert_block(Matrix(1, 1, 1), C, D);
            create_down.insert_block(Matrix(1, 1, 1), A, C);
            create_down.insert_block(Matrix(1, 1, 1), B, D);
            
            destroy_up.insert_block(Matrix(1, 1, 1), B, A);
            destroy_up.insert_block(Matrix(1, 1, 1), D, C);
            destroy_down.insert_block(Matrix(1, 1, 1), C, A);
            destroy_down.insert_block(Matrix(1, 1, 1), D, B);
            
            count_up.insert_block(Matrix(1, 1, 1), B, B);
            count_up.insert_block(Matrix(1, 1, 1), D, D);
            count_down.insert_block(Matrix(1, 1, 1), C, C);
            count_down.insert_block(Matrix(1, 1, 1), D, D);
            
            doubly_occ.insert_block(Matrix(1, 1, 1), D, D);
            
            sign_up.insert_block(Matrix(1, 1, 1), A, A);
            sign_up.insert_block(Matrix(1, 1, -1), B, B);
            sign_up.insert_block(Matrix(1, 1, 1), C, C);
            sign_up.insert_block(Matrix(1, 1, -1), D, D);
            
            sign_down.insert_block(Matrix(1, 1, 1), A, A);
            sign_down.insert_block(Matrix(1, 1, 1), B, B);
            sign_down.insert_block(Matrix(1, 1, -1), C, C);
            sign_down.insert_block(Matrix(1, 1, -1), D, D);
            
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
                        term.fill_operator = sign_up;
                        gemm(create_up, sign_up, tmp);
                        term.operators.push_back( std::make_pair(p, -ti*tmp) );
                        term.operators.push_back( std::make_pair(*hopto, destroy_up) );
                        terms.push_back(term);
                    }
                    { // t*c_up*cdag_up
                        hamterm_t term;
                        term.fill_operator = sign_up;
                        gemm(sign_up, destroy_up, tmp);
                        term.operators.push_back( std::make_pair(p, -ti*tmp) );
                        term.operators.push_back( std::make_pair(*hopto, create_up) );
                        terms.push_back(term);
                    }
                    { // t*cdag_down*c_down
                        hamterm_t term;
                        term.fill_operator = sign_down;
                        gemm(create_down, sign_down, tmp);
                        term.operators.push_back( std::make_pair(p, -ti*tmp) );
                        term.operators.push_back( std::make_pair(*hopto, destroy_down) );
                        terms.push_back(term);
                    }
                    { // t*c_down*cdag_down
                        hamterm_t term;
                        term.fill_operator = sign_down;
                        gemm(sign_down, destroy_down, tmp);
                        term.operators.push_back( std::make_pair(p, -ti*tmp) );
                        term.operators.push_back( std::make_pair(*hopto, create_down) );
                        terms.push_back(term);
                    }
                }
            }
        }
        
        int n_terms() const
        {
            return terms.size();
        }
        hamterm_t operator[](int i) const
        {
            return terms[i];
        }
        
        op_t get_identity() const
        {
            return ident;
        }
        
        Index<TwoU1> get_phys() const
        {
            Index<TwoU1> phys;
            
            TwoU1::charge A(0), B(0), C(0), D(1);
            B[0]=1; C[1]=1;
            
            phys.insert(std::make_pair(A, 1));
            phys.insert(std::make_pair(B, 1));
            phys.insert(std::make_pair(C, 1));
            phys.insert(std::make_pair(D, 1));
            
            return phys;
        }
        
    private:
        std::vector<hamterm_t> terms;
        op_t ident;
        
        double get_t (BaseParameters & parms, int i)
        {
            std::ostringstream key;
            key << "t" << (i+1);
            return (parms.is_set(key.str())) ? parms.get<double>(key.str()) : parms.get<double>("t");
        }
    };
    
    
} // namespace

#endif
