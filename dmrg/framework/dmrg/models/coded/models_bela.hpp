/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef BELA_HAMILTONIANS_H
#define BELA_HAMILTONIANS_H

#include <sstream>

#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"

namespace app {
    
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
            
            std::cout << "Wrong guy." << std::endl;
            exit(0);
        }
        
        void do_init(const Lattice & lat, BaseParameters & parms, std::complex<double>)
        {
            using std::cout;
            using std::endl;
            
            phys.insert(std::make_pair(-1, 1));
            phys.insert(std::make_pair(1, 1));
            
            ident.insert_block(Matrix(1, 1, 1), -1, -1);
            ident.insert_block(Matrix(1, 1, 1), 1, 1);
            
            splus.insert_block(Matrix(1, 1, 1), -1, 1);
            sminus.insert_block(Matrix(1, 1, 1), 1, -1);
            
            sz.insert_block(Matrix(1, 1, 0.5), 1, 1);
            sz.insert_block(Matrix(1, 1, -0.5), -1, -1);
            
            double
            J1 = parms.get<double>("J1")
            , J2 = parms.get<double>("J2")
            , delta = parms.get<double>("delta")
            , KK = parms.get<double>("K2")
            , funny = parms.get<int>("funny");
            
            int L = lat.size();
            std::vector<int> last_adj = lat.all(L-1);
            bool open = !(std::count(last_adj.begin(), last_adj.end(), 0) > 0);
            if (open)
                cout << "OBC" << endl;
            else
                cout << "PBC" << endl;

#define NN(p,op1,q,op2) \
{ hamterm_t term; \
term.fill_operator = ident; \
term.operators.push_back( make_pair(p, op1) ); \
term.operators.push_back( make_pair(q, op2) ); \
terms.push_back(term); }

#define GRP(p1,p2,J) \
NN(p1, J/2*splus, p2, sminus); \
NN(p1, J/2*sminus, p2, splus); \
NN(p1, delta*J*sz, p2, sz);

            if (KK == 0) {

                // Heisenberg term
                for (int p = 0; p < L; ++p) {
                    if (open && p+1 >= L)
                        continue;
                    NN(p, J1/2*splus, (p+1)%L, sminus);
                    NN(p, J1/2*sminus, (p+1)%L, splus);
                    NN(p, delta*J1*sz, (p+1)%L, sz);
                }

                // NNN HB term
                for (int p = 0; p < L - 2; ++p) {
                    if (open && p+2 >= L)
                        continue;
                    NN(p, J2/2*splus, (p+2)%L, sminus);
                    NN(p, J2/2*sminus, (p+2)%L, splus);
                    NN(p, delta*J2*sz, (p+2)%L, sz);
                }

            } else {

                for (int p = 0; p < L; ++p) {
                    GRP(p, (p+3)%L, J2);
                    if (p % 3 == 2)
                        continue;
                    GRP(p, (p+1)%L, J1);
                }

                for (int p = 1; p < L; p += 3) {
                    GRP(p, (p+2)%L, J1);
                    GRP(p, (p+4)%L, J1);
                }
            }
#undef NN
#undef GRP

#define NNN(p,op1,op2,op3) \
{ hamterm_t term; \
term.fill_operator = ident; \
term.operators.push_back( make_pair(p, op1) ); \
term.operators.push_back( make_pair((p+1)%L, op2) ); \
term.operators.push_back( make_pair((p+2)%L, op3) ); \
terms.push_back(term); }

            if (KK == 0) {

                for (int p = 0; p < L; ++p) {
                    if (open && p+2 >= L)
                        continue;

                    double K = (p % 2 == 0 ? parms.get<double>("K0") : parms.get<double>("K1"));
                    if (funny == 1) {
                        int p1 = L/2-1;
                        if (p == p1)
                            continue;
                    } else if (funny == 2) {
                        int p1 = L/2-1;
                        K = parms.get<double>("K0");
                        if (p < p1 && p1 % 2 == 1)
                            K *= -1;
                    } else if (funny == 3) {
                        int p1 = L/2-1;
                        K = parms.get<double>("K0") * powf(-1, p);
                        if (p == p1 || p == p1-1 || p == p1+1)
                            K = parms.get<double>("K1") * powf(-1, p);
                        //else if (p == p1 - 1 || p == p1 + 1)
                        //    K = parms.get<double>("K1") * powf(-1, p);
                    } else if (funny == 3) {
                        int p1 = L/2-1;
                        K = parms.get<double>("K0") * powf(-1, p);
                        if (p == p1)
                            continue;
                        else if (p == p1 - 1 || p == p1 + 1)
                            K = parms.get<double>("K1") * powf(-1, p);
                    }

                    std::complex<double> f(0, 0.5);

                    NNN(p,  f*K * splus, sminus, delta*sz);
                    NNN(p, -f*K * splus, sz, sminus);

                    NNN(p,  f*K * sminus, sz, splus);
                    NNN(p, -f*K * sminus, splus, delta*sz);

                    NNN(p,  f*K * delta*sz, splus, sminus);
                    NNN(p, -f*K * delta*sz, sminus, splus);
                }

            } else

#undef NNN

#define NNN(p1,op1,p2,op2,p3,op3) \
{ hamterm_t term; \
term.fill_operator = ident; \
term.operators.push_back( make_pair(p1, op1) ); \
term.operators.push_back( make_pair(p2, op2) ); \
term.operators.push_back( make_pair(p3, op3) ); \
terms.push_back(term); }

#define TERM(p1_,p2_,p3_) \
if (true) { \
cout << "Operator on " << p1_ << " " << p2_ << " " << p3_ << endl; \
int p1=(p1_)%L, p2=(p2_)%L, p3=(p3_)%L; \
NNN(p1,  f*KK * splus, p2, sminus, p3, delta*sz); \
NNN(p1, -f*KK * splus, p2, sz, p3, sminus); \
NNN(p1,  f*KK * sminus, p2, sz, p3, splus); \
NNN(p1, -f*KK * sminus, p2, splus, p3, delta*sz); \
NNN(p1,  f*KK * delta*sz, p2, splus, p3, sminus); \
NNN(p1, -f*KK * delta*sz, p2, sminus, p3, splus); }

            { // see above: KK != 0

                for (int p = 0; p < L; ++p) {
                    std::complex<double> f(0, 0.5);
                    /*if (p % 3 == 0) {
                        TERM(p, p+1, p+3);
                    } else if (p % 3 == 1) {
                        TERM(p, p+2, p+3);
                        TERM(p, p+4, p+3);
                    } else if (p % 3 == 2) {
                        TERM(p, p-1, p+3);
                    }*/
                    if (p % 3 == 0) {
                        TERM(p, p+1, p+3);
                    } else if (p % 3 == 1) {
                        TERM(p, p+2, p+3);
                        TERM(p, p+1, p+3);
                    } else if (p % 3 == 2) {
                        TERM(p, p+2, p+3);
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

        Measurements<Matrix, U1> measurements () const
        {
            return Measurements<Matrix, U1>();
        }
        
    private:
        op_t ident, splus, sminus, sz;
        Index<U1> phys;
        
        std::vector<hamterm_t> terms;
    };
        
} // namespace

#endif
