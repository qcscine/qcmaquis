#ifndef B_2U1_HAMILTONIANS_H
#define B_2U1_HAMILTONIANS_H

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "b_generate_mpo.h"

#include <vector>
#include <utility>

namespace b_mpos {
    template<class Matrix>
    class TwoU1_Spin1BlBq : public Hamiltonian<Matrix, TwoU1>
    {
    public:
        typedef typename Hamiltonian<Matrix, TwoU1>::op_t op_t;
        typedef std::pair<op_t, op_t> op_pair;
        typedef std::vector<op_pair> op_pairs;
        
        TwoU1_Spin1BlBq(BaseParameters & parms_)
        : Jbl(cos(M_PI * parms_.get<double>("theta")))
        , Jbq(sin(M_PI * parms_.get<double>("theta")))
        , parms(parms_)
        {  
            block_matrix<Matrix, TwoU1> ident, aa, bb, cc, ab, ba, ac, ca, bc, cb;
            
            TwoU1::charge A, B, C;
            B[0] = 1;
            C[1] = 1;
            
#define define_op(name, I, J) name.insert_block(Matrix(1, 1, 1), I, J)
            
            define_op(ident, A, A);
            define_op(ident, B, B);
            define_op(ident, C, C);
            
            define_op(aa, A, A);
            define_op(bb, B, B);
            define_op(cc, C, C);
            
            define_op(ab, A, B);
            define_op(ba, B, A);
            
            define_op(ac, A, C);
            define_op(ca, C, A);
            
            define_op(bc, B, C);
            define_op(cb, C, B);
            
#undef define_op
            
#define term(a,b) ops.push_back(make_pair(a, b))
            
            term(ident, ident);
            
            term(aa, aa);
            term(bb, bb);
            term(cc, cc);
            
            term(ab, ba);
            term(ac, ca);
            term(bc, cb);
            
            term(ba, ab);
            term(ca, ac);
            term(cb, bc);
            
#undef term
        }
        
        op_t get_identity() { return ops[0].first; }
        op_t get_free() { return ops[0].first; }
        
        int num_2site_ops() { return 0; }
        op_pair get_2site_op(int i) { return op_pair(); }
        
        int num_1site_ops() { return 0; }
        op_t get_1site_op(int) { return op_t(); }
        
        Index<TwoU1> get_phys()
        {
            Index<TwoU1> phys;
            
            TwoU1::charge A, B, C;
            B[0] = 1;
            C[1] = 1;
            
            phys.insert(std::make_pair(A, 1));
            phys.insert(std::make_pair(B, 1));
            phys.insert(std::make_pair(C, 1));
            return phys;
        }
        
        void push_extra_terms(MPOMaker<Matrix, TwoU1> & mm, b_adj::Adjacency& adj)
        {
            std::vector<std::pair<std::size_t, op_t> > terms;
            
            for (int p = 0; p < adj.size(); ++p) {
                vector<int> neighs = adj.forward(p);
                for (vector<int>::iterator neigh = neighs.begin(); neigh != neighs.end(); ++neigh) {
                    double K;
                    if (adj.bond_type(p, *neigh) == 0)
                        K = parms.get<double>("K0");
                    else
                        K = parms.get<double>("K1");
                    
                    if (K == 0)
                        continue;
                    
                    maquis::cout << p << " " << *neigh << " K=" << K << std::endl;
                    
                    for (int i = 0; i < ops.size()-1; ++i)
                    {
                        terms.clear();
                        terms.push_back( make_pair( p, K*ops[i+1].first ) );
                        terms.push_back( make_pair( *neigh, ops[i+1].second ) );
                        mm.add_term(terms);
                    }
                }
            }
            
            block_matrix<Matrix, TwoU1> c[3];
            
            TwoU1::charge A, B, C;
            B[0] = 1;
            C[1] = 1;
            
            c[0].insert_block(Matrix(1, 1, 1), A, A);
            c[1].insert_block(Matrix(1, 1, 1), B, B);
            c[2].insert_block(Matrix(1, 1, 1), C, C); 

            double h0 = parms.get<double>("h0");
            int W = parms.get<int>("W");
            
            if (parms.get<int>("pin") > 0) {
                for (int i = 0; i < parms.get<int>("pin"); ++i)
                {
                    terms.clear();
                    terms.push_back( std::make_pair(i, h0*c[(i+2*(i/W))%3]) );
                    mm.add_term(terms);
                }
            }
            
            if (parms.get<int>("pin") < 0) {
                for (int i = 0; i < adj.size(); ++i) {
                    if (adj.site_type(i) > 0) {
                        maquis::cout << "Pinning on site " << i << " to color " << adj.site_type(i)-1 << std::endl;
                        
                        terms.clear();
                        terms.push_back(std::make_pair(i, h0*c[adj.site_type(i)-1]));
                        mm.add_term(terms);
                    }
                }
            }
        }
        
    private:
        double Jbl, Jbq;
        op_pairs ops;
        BaseParameters & parms;
    };
    
    template<class Matrix>
    class TwoU1_FermiHubbard : public Hamiltonian<Matrix, TwoU1>
    {
        typedef typename Hamiltonian<Matrix, TwoU1>::op_t op_t;
        typedef std::pair<op_t, op_t> op_pair;
        typedef std::vector<op_pair> op_pairs;
        
    public:
        TwoU1_FermiHubbard(double t_, double U_) : t(t_), U(U_)
        {
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
        }
        
        void push_extra_terms(MPOMaker<Matrix, TwoU1> & mm, b_adj::Adjacency& adj)
        {
            vector<pair<size_t, op_t> > term;
            
            // U term
            for (int p = 0; p < adj.size(); ++p)
            {
                term.push_back(make_pair(p, U*doubly_occ));
                mm.add_term(term);
                term.clear();
            }
            
            block_matrix<Matrix, TwoU1> tmp;
            
            for (int p = 0; p < adj.size(); ++p)
            {
                vector<int> neighs = adj.forward(p);
                
                for (vector<int>::iterator hopto = neighs.begin();
                     hopto != neighs.end(); ++hopto)
                {
                    term.clear();
                    
                    for (int f0 = std::min(p, *hopto)+1;
                         f0 < std::max(p, *hopto); ++f0)
                        term.push_back(make_pair(f0, sign_up));
                    
                    gemm(create_up, sign_up, tmp);
                    term.push_back(make_pair(p, -t*tmp));
                    term.push_back(make_pair(*hopto, destroy_up));
                    
                    mm.add_term(term);
                    term.clear();
                    
                    for (int f0 = std::min(p, *hopto)+1;
                         f0 < std::max(p, *hopto); ++f0)
                        term.push_back(make_pair(f0, sign_up));
                    
                    gemm(sign_up, destroy_up, tmp);
                    term.push_back(make_pair(p, -t*tmp));
                    term.push_back(make_pair(*hopto, create_up));
                    
                    mm.add_term(term);
                    term.clear();
                    
                    for (int f0 = std::min(p, *hopto)+1;
                         f0 < std::max(p, *hopto); ++f0)
                        term.push_back(make_pair(f0, sign_down));
                    
                    gemm(create_down, sign_down, tmp);
                    term.push_back(make_pair(p, -t*tmp));
                    term.push_back(make_pair(*hopto, destroy_down));
                    
                    mm.add_term(term);
                    term.clear();
                    
                    for (int f0 = std::min(p, *hopto)+1;
                         f0 < std::max(p, *hopto); ++f0)
                        term.push_back(make_pair(f0, sign_down));
                    
                    gemm(sign_down, destroy_down, tmp);
                    term.push_back(make_pair(p, -t*tmp));
                    term.push_back(make_pair(*hopto, create_down));
                    
                    mm.add_term(term);
                    term.clear();
                }
            }
        }
        
        op_t get_identity() { return ident; }
        op_t get_free() { return ident; }
        
        int num_2site_ops() { return 0; }
        op_pair get_2site_op(int) { return op_pair(); }
        
        int num_1site_ops() { return 0; }
        op_t get_1site_op(int) { return op_t(); }
        
        Index<TwoU1> get_phys()
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
        double t, U;
        op_pairs ops;
        
        block_matrix<Matrix, TwoU1> create_up, create_down;
        block_matrix<Matrix, TwoU1> destroy_up, destroy_down;
        block_matrix<Matrix, TwoU1> count_up, count_down;            
        block_matrix<Matrix, TwoU1> ident, doubly_occ;
        block_matrix<Matrix, TwoU1> sign_up, sign_down;
    };
}

#endif
