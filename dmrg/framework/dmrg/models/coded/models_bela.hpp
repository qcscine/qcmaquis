/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MODELS_CODED_BELA_H
#define MODELS_CODED_BELA_H

#include <sstream>

#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"

/* Chiral spins */
template<class Matrix>
class Chiral : public Model<Matrix, U1>
{
    typedef Hamiltonian<Matrix, U1> ham;        
    typedef typename ham::hamterm_t hamterm_t;        
    typedef typename ham::op_t op_t;
    typedef Measurement_Term<Matrix, U1> mterm_t;
    
public:
    // This can only work for complex types, that's why there's a callback to do_init with a tag
    Chiral(const Lattice & lat, BaseParameters & parms)
    {
        do_init(lat, parms, typename Matrix::value_type());
    }
    
    void do_init(const Lattice & lat, BaseParameters & parms, double)
    {
        phys.insert(std::make_pair(0, 1));
        phys.insert(std::make_pair(1, 1));
        
        maquis::cout << "Wrong guy." << std::endl;
        exit(0);
    }
    
    void do_init(const Lattice & lat, BaseParameters & parms, std::complex<double>)
    {
        phys.insert(std::make_pair(-1, 1));
        phys.insert(std::make_pair(1, 1));
        
        ident.insert_block(Matrix(1, 1, 1), -1, -1);
        ident.insert_block(Matrix(1, 1, 1), 1, 1);
        
        splus.insert_block(Matrix(1, 1, 1), -1, 1);
        sminus.insert_block(Matrix(1, 1, 1), 1, -1);
        
        sz.insert_block(Matrix(1, 1, 0.5), 1, 1);
        sz.insert_block(Matrix(1, 1, -0.5), -1, -1);
        
        K2 = parms.get<double>("K2");
        K3 = parms.get<double>("K3");
        K4 = parms.get<double>("K4");
        delta = 1;
        stag = parms.get<int>("stag");
        alt = parms.get<int>("alt");
        
        L = lat.size();
        std::vector<int> last_adj = lat.all(L-1);
        open = !(std::count(last_adj.begin(), last_adj.end(), 0) > 0);
        if (open)
            maquis::cout << "OBC" << std::endl;
        else
            maquis::cout << "PBC" << std::endl;
        if (stag)
            maquis::cout << "Staggered" << std::endl;
        else
            maquis::cout << "Non-staggered" << std::endl;
        
        if (K3 == 0 && K4 == 0 && K2 != 0) {
            f = std::complex<double>(0, K2 * 0.5);
            for (int p = 0; p < L; ++p) {
                if (stag) {
                    if (p % 3 == 0) {
                        TERM(p, p+3, p+1);
                    } else if (p % 3 == 1) {
                        TERM(p, p+1, p+4);
                    }
                } else {
                    if (p % 6 == 0) {
                        TERM(p, p+3, p+1);
                    } else if (p % 6 == 1) {
                        TERM(p, p+4, p+1);
                    } else if (p % 6 == 3) {
                        TERM(p, p+1, p+3);
                    } else if (p % 6 == 4) {
                        TERM(p, p+1, p+4);
                    }
                }
                
                if (open && p == L - 4) {
                    TERM(p, p+1, p+3);
                }
            }
        } else if (K3 != 0 && K4 == 0 && K2 == 0) {
            if (alt == 0) {
                std::cout << "Setting up operators for K3 non-alternative." << std::endl;
                f = std::complex<double>(0, K3 * 0.5);
                int pos[] = {
                    0, 2, 1,
                    2, 3, 4,
                    4, 5, 15,
                    5, 7, 18,
                    6, 8, 7,
                    8, 9, 10,
                    10, 11, 21,
                    11, 1, 12
                };
                std::cout << L << " " << L/12 << std::endl;
                int lmax = int(L/12)*12;
                for (int p = 0; p < L; p += 12) {
                    std::cout << p << std::endl;
                    for (int k = 0; k < 8; ++k) {
                        if (p+pos[k*3+0] < lmax && p+pos[k*3+1] < lmax && p+pos[k*3+2] < lmax) {
                            if (stag && k % 2 == 1)
                                TERM(p+pos[k*3+0], p+pos[k*3+2], p+pos[k*3+1]);
                            else
                                TERM(p+pos[k*3+0], p+pos[k*3+1], p+pos[k*3+2]);
                        }
                    }
                }
            } else {
                std::cout << "Setting up operators for K3 alternative." << std::endl;
                f = std::complex<double>(0, K3 * 0.5);
                int pos[] = {
                    4, 3, 2,
                    0, 2, 16,
                    1, 0, 17,
                    9, 5, 1,
                    16, 15, 14,
                    12, 14, 4,
                    13, 12, 5,
                    21, 17, 13,
                    10, 9, 8,
                    6, 8, 22,
                    7, 6, 23,
                    27, 11, 7,
                    22, 21, 20,
                    18, 20, 10,
                    19, 18, 11,
                    39, 23, 19
                };
                
                std::cout << L << " " << L/24 << std::endl;
                for (int p = 0; p < L; p += 24) {
                    std::cout << p << std::endl;

                    for (int k = 0; k < 16; ++k) {
                        TERM(p+pos[k*3+0], p+pos[k*3+1], p+pos[k*3+2]);
                    }
                }
                
                int last = (L/24-1)*24;
                TERM(L-2, last+11, last+7);
                TERM(L-1, last+23, last+19);
            }
            if (open) {
                int p0 = (int(L/12)-1)*12;
                TERM(p0+4, p0+5, p0+13);
                TERM(p0+10, p0+11, p0+15);
                if (stag) {
                    TERM(p0+7, p0+5, p0+14);
                    TERM(p0+1, p0+11, p0+12);
                } else {
                    TERM(p0+5, p0+7, p0+14);
                    TERM(p0+1, p0+12, p0+11);
                }
            }
        } else if (K4 != 0) {
            f = std::complex<double>(0, K4 * 0.5);
            int pos_nonstag[] = {
                0, 3, 1,
                1, 4, 2,
                3, 5, 7,
                4, 9, 6,
            };
            int pos_stag[] = {
                0, 3, 1,
                1, 2, 4,
                3, 7, 5,
                4, 9, 6
            };
            
            int * pos;
            if (stag)
                pos = &pos_stag[0];
            else
                pos = &pos_nonstag[0];
            
            for (int p = 0; p < L; p += 7) {
                std::cout << p << std::endl;
                for (int k = 0; k < 4; ++k) {
                    TERM(p+pos[k*3+0], p+pos[k*3+1], p+pos[k*3+2]);
                }
            }
        }
        

#undef TERM
#undef NNN
    }
    
    Hamiltonian<Matrix, U1> H () const
    {        
        return ham(phys, ident, terms);
    }
    
    Index<U1> get_phys() const
    {
        return phys;
    }
    
    Measurements<Matrix, U1> measurements () const
    {
        // return Measurements<Matrix, U1>();
        Measurements<Matrix, U1> ret;
        for (typename std::vector<boost::shared_ptr<mterm_t> >::const_iterator it = m_terms.begin(); it != m_terms.end(); ++it)
            ret.add_term(**it);
        ret.set_identity(ident);
        return ret;
    }
    
private:
    void TERM(int p1__, int p2__, int p3__)
    {
#define NNN(p1,op1,p2,op2,p3,op3) \
{ hamterm_t term; \
term.fill_operator = ident; \
term.operators.push_back( make_pair(p1, op1) ); \
term.operators.push_back( make_pair(p2, op2) ); \
term.operators.push_back( make_pair(p3, op3) ); \
terms.push_back(term); }

#define NNN_m(p1,op1,p2,op2,p3,op3) \
{ std::vector< std::pair<int, op_t> > t__; \
t__.push_back( make_pair(p1, op1) ); \
t__.push_back( make_pair(p2, op2) ); \
t__.push_back( make_pair(p3, op3) ); \
measterm->custom_ops.push_back(t__); }

#define TERM__(p1_,p2_,p3_) \
if (!open || (p1_ < L && p2_ < L && p3_ < L)) { \
maquis::cout << "Operator on " << p1_ << " " << p2_ << " " << p3_ << std::endl; \
int p1=(p1_)%L, p2=(p2_)%L, p3=(p3_)%L; \
{ NNN(p1,  f * splus, p2, sminus, p3, delta*sz); \
NNN(p1, -f * splus, p2, sz, p3, sminus); \
NNN(p1,  f * sminus, p2, sz, p3, splus); \
NNN(p1, -f * sminus, p2, splus, p3, delta*sz); \
NNN(p1,  f * delta*sz, p2, splus, p3, sminus); \
NNN(p1, -f * delta*sz, p2, sminus, p3, splus); } \
{ boost::shared_ptr<mterm_t> measterm(new mterm_t()); \
measterm->fill_operator = ident; \
std::ostringstream name_oss; \
name_oss << "BE " << p1 << " " << p2 << " " << p3; \
measterm->name = name_oss.str(); \
measterm->type = mterm_t::Custom; \
NNN_m(p1,  f * splus, p2, sminus, p3, delta*sz); \
NNN_m(p1, -f * splus, p2, sz, p3, sminus); \
NNN_m(p1,  f * sminus, p2, sz, p3, splus); \
NNN_m(p1, -f * sminus, p2, splus, p3, delta*sz); \
NNN_m(p1,  f * delta*sz, p2, splus, p3, sminus); \
NNN_m(p1, -f * delta*sz, p2, sminus, p3, splus); \
m_terms.push_back(boost::shared_ptr<mterm_t>(measterm)); } }

        TERM__(p1__, p2__, p3__);
    }
    
    op_t ident, splus, sminus, sz;
    Index<U1> phys;
    
    std::vector<hamterm_t> terms;
    std::vector<boost::shared_ptr<mterm_t> > m_terms;
    
    double K2, K3, K4, delta;
    int stag, L, alt;
    bool open;
    std::complex<double> f;
};

#endif
