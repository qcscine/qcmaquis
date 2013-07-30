#ifndef B_HAMILTONIANS_H
#define B_HAMILTONIANS_H

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "b_generate_mpo.h"
#include "b_DmrgParameters.h"

#include <vector>
#include <utility>

namespace b_mpos {
    template<class Matrix>
    class Heisenberg : public Hamiltonian<Matrix, U1>
    {
    public:
        typedef typename Hamiltonian<Matrix, U1>::op_t op_t;
        typedef std::pair<op_t, op_t> op_pair;
        typedef std::vector<op_pair> op_pairs;
        
        Heisenberg(double Jxy_, double Jz_) : Jxy(Jxy_), Jz(Jz_)
        {
            block_matrix<Matrix, U1> ident, splus, sminus, sz;
            
            ident.insert_block(Matrix(1, 1, 1), -1, -1);
            ident.insert_block(Matrix(1, 1, 1), 1, 1);
            
            splus.insert_block(Matrix(1, 1, 1), -1, 1);
            
            sminus.insert_block(Matrix(1, 1, 1), 1, -1);
            
            sz.insert_block(Matrix(1, 1, 0.5), 1, 1);
            sz.insert_block(Matrix(1, 1, -0.5), -1, -1);
            
            ops.push_back(make_pair(ident, ident));
            ops.push_back(make_pair(Jxy/2*splus, sminus));
            ops.push_back(make_pair(Jxy/2*sminus, splus));
            ops.push_back(make_pair(Jz*sz, sz));
        }
        
        op_t get_identity() { return ops[0].first; }
        op_t get_free() { return ops[0].first; }
        
        int num_2site_ops() { return ops.size()-1; }
        op_pair get_2site_op(int i) { return ops[i+1]; }
        
        int num_1site_ops() { return 0; }
        op_t get_1site_op(int) { return op_t(); }
        
        Index<U1> get_phys()
        {
            Index<U1> phys;
            phys.insert(std::make_pair(1, 1));
            phys.insert(std::make_pair(-1, 1));
            return phys;
        }
        
    private:
        double Jxy, Jz;
        op_pairs ops;
    };
    
    template<class Matrix>
    class HCB : public Hamiltonian<Matrix, U1>
    {
        typedef typename Hamiltonian<Matrix, U1>::op_t op_t;
        typedef std::pair<op_t, op_t> op_pair;
        typedef std::vector<op_pair> op_pairs;
        
    public:
        HCB()
        {   
            block_matrix<Matrix, U1> create, destroy, ident, count;
            
            ident.insert_block(Matrix(1, 1, 1), 0, 0);
            ident.insert_block(Matrix(1, 1, 1), 1, 1);
            
            create.insert_block(Matrix(1, 1, 1), 0, 1);
            destroy.insert_block(Matrix(1, 1, 1), 1, 0);
            
            count.insert_block(Matrix(1, 1, 1), 1, 1);
            
            double t = 1;
            
#define term(a,b) ops.push_back(make_pair(a, b))
            term(ident, ident);
            term(-t*create, destroy);
            term(-t*destroy, create);
#undef term
        }
        
        op_t get_identity() { return ops[0].first; }
        op_t get_free() { return ops[0].first; }
        
        int num_2site_ops() { return ops.size()-1; }
        op_pair get_2site_op(int i) { return ops[i+1]; }
        
        int num_1site_ops() { return 0; }
        op_t get_1site_op(int) { return op_t(); }
        
        Index<U1> get_phys()
        {
            Index<U1> phys;
            phys.insert(std::make_pair(0, 1));
            phys.insert(std::make_pair(1, 1));
            return phys;
        }
        
    private:
        op_pairs ops;
    };
    
    template<class Matrix>
    class FreeFermions : public Hamiltonian<Matrix, U1>
    {
        typedef typename Hamiltonian<Matrix, U1>::op_t op_t;
        typedef std::pair<op_t, op_t> op_pair;
        typedef std::vector<op_pair> op_pairs;
        
    public:
        FreeFermions()
        {   
            block_matrix<Matrix, U1> create, destroy;
            create.insert_block(Matrix(1, 1, 1), 0, 1);
            destroy.insert_block(Matrix(1, 1, 1), 1, 0);
            
            double t = 1;
            
#define term(a,b) ops.push_back(make_pair(a, b))
            term(-t*create, destroy);
            term(-t*destroy, create);
#undef term
        }
        
        op_t get_identity()
        { 
            block_matrix<Matrix, U1> ident;
            
            ident.insert_block(Matrix(1, 1, 1), 0, 0);
            ident.insert_block(Matrix(1, 1, 1), 1, 1);
            
            return ident;
        }
        
        op_t get_free()
        {
            block_matrix<Matrix, U1> sign;
            
            sign.insert_block(Matrix(1, 1, 1), 0, 0);
            sign.insert_block(Matrix(1, 1, -1), 1, 1);
            
            return sign;
        }
        
        int num_2site_ops() { return ops.size(); }
        op_pair get_2site_op(int i) { return ops[i]; }
        
        int num_1site_ops() { return 0; }
        op_t get_1site_op(int) { return op_t(); }
        
        Index<U1> get_phys()
        {
            Index<U1> phys;
            phys.insert(std::make_pair(0, 1));
            phys.insert(std::make_pair(1, 1));
            return phys;
        }
        
    private:
        op_pairs ops;
    };
    
    template<class Matrix>
    class AlternateFreeFermions : public Hamiltonian<Matrix, U1>
    {
        typedef typename Hamiltonian<Matrix, U1>::op_t op_t;
        typedef std::pair<op_t, op_t> op_pair;
        typedef std::vector<op_pair> op_pairs;
        
    public:
        AlternateFreeFermions(BaseParameters & parms_)
        : parms(parms_)
        {
            create.insert_block(Matrix(1, 1, 1), 0, 1);
            destroy.insert_block(Matrix(1, 1, 1), 1, 0);
            
            count.insert_block(Matrix(1, 1, 1), 1, 1);
            omcount.insert_block(Matrix(1, 1, 1), 0, 0);
            
            ident.insert_block(Matrix(1, 1, 1), 0, 0);
            ident.insert_block(Matrix(1, 1, 1), 1, 1);
            
            sign.insert_block(Matrix(1, 1, 1), 0, 0);
            sign.insert_block(Matrix(1, 1, -1), 1, 1);
        }
        
        op_t get_identity() { return ident; }
        op_t get_free() { return ident; }
        
        int num_2site_ops() { return 0; }
        op_pair get_2site_op(int) { return op_pair(); }
        int num_1site_ops() { return 0; }
        op_t get_1site_op(int) { return op_t(); }
        
        void push_extra_terms(MPOMaker<Matrix, U1> & mm,
                              b_adj::Adjacency & adj)
        {
            for (int p = 0; p < adj.size(); ++p)
            {   
                vector<pair<size_t, op_t> > term;
                
                vector<int> neighs = adj.all(p);
                
                // hopping
                for (vector<int>::iterator hopto = neighs.begin();
                     hopto != neighs.end(); ++hopto)
                {
                    double phase = 1;
                    if (parms["twist"] == 1 && adj.wraps_pbc(p, *hopto))
                        phase = -1;
                    
                    term.clear();
                    
                    for (int f0 = std::min(p, *hopto)+1;
                         f0 < std::max(p, *hopto); ++f0)
                        term.push_back(make_pair(f0, sign));
                    
                    term.push_back(make_pair(p, -phase*destroy));
                    term.push_back(make_pair(*hopto, create));
                    
                    mm.add_term(term);
                }
            }
        }
        
        Index<U1> get_phys()
        {
            Index<U1> phys;
            phys.insert(std::make_pair(0, 1));
            phys.insert(std::make_pair(1, 1));
            return phys;
        }
        
    private:
        int sweep;
        BaseParameters & parms;
        
        block_matrix<Matrix, U1> create, destroy, count, omcount, ident, sign;
    };
    
    template<class Matrix>
    class Spin1BlBq : public Hamiltonian<Matrix, U1>
    {
    public:
        typedef typename Hamiltonian<Matrix, U1>::op_t op_t;
        typedef std::pair<op_t, op_t> op_pair;
        typedef std::vector<op_pair> op_pairs;
        
        Spin1BlBq(BaseParameters & parms_)
        : Jbl(cos(M_PI * parms_["theta"]))
        , Jbq(sin(M_PI * parms_["theta"]))
        , parms(parms_)
        {   
            block_matrix<Matrix, U1> ident, splus, sminus, sz, spp, smm, spm, smp, szz, szp, spz, szm, smz;
            
            ident.insert_block(Matrix(1, 1, 1), -1, -1);
            ident.insert_block(Matrix(1, 1, 1), 0, 0);
            ident.insert_block(Matrix(1, 1, 1), 1, 1);
            
            splus.insert_block(Matrix(1, 1, 1), -1, 0);
            splus.insert_block(Matrix(1, 1, 1), 0, 1);
            
            sminus.insert_block(Matrix(1, 1, 1), 1, 0);
            sminus.insert_block(Matrix(1, 1, 1), 0, -1);
            
            sz.insert_block(Matrix(1, 1, 1), 1, 1);
            sz.insert_block(Matrix(1, 1, 0), 0, 0);
            sz.insert_block(Matrix(1, 1, -1), -1, -1);
            
            gemm(splus, splus, spp);
            gemm(sminus, sminus, smm);
            gemm(splus, sminus, spm);
            gemm(sminus, splus, smp);
            gemm(sz, sz, szz);
            gemm(sz, splus, szp);
            gemm(splus, sz, spz);
            gemm(sz, sminus, szm);
            gemm(sminus, sz, smz);
            
#define term(a,b) ops.push_back(make_pair(a, b))
            
            term(ident, ident);
            term(Jbl*splus, sminus);
            term(Jbl*sminus, splus);
            term(Jbl*sz, sz);
            
            term(Jbq*spp, smm);
            term(Jbq*smm, spp);
            term(Jbq*szz, szz);
            
            term(Jbq*spm, smp);
            term(Jbq*smp, spm);
            
            term(Jbq*spz, smz);
            term(Jbq*szp, szm);
            
            term(Jbq*smz, spz);
            term(Jbq*szm, szp);
            
#undef term
        }
        
        op_t get_identity() { return ops[0].first; }
        op_t get_free() { return ops[0].first; }
        
        int num_2site_ops() { return ops.size()-1; }
        op_pair get_2site_op(int i) { return ops[i+1]; }
        
        int num_1site_ops() { return 0; }
        op_t get_1site_op(int) { return op_t(); }
        
        Index<U1> get_phys()
        {
            Index<U1> phys;
            phys.insert(std::make_pair(1, 1));
            phys.insert(std::make_pair(0, 1));
            phys.insert(std::make_pair(-1, 1));
            return phys;
        }
        
        void push_extra_terms(MPOMaker<Matrix, U1> & mm, b_adj::Adjacency& adj)
        {
            maquis::cout << "Adding extra term." << std::endl;
            std::vector<std::pair<std::size_t, op_t> > terms;
            block_matrix<Matrix, U1> c[3];
            
            c[0].insert_block(Matrix(1, 1, 1), -1, -1);
            c[1].insert_block(Matrix(1, 1, 1), 0, 0);
            c[2].insert_block(Matrix(1, 1, 1), 1, 1);
            
            double h0 = parms["h0"];
            int W = parms["W"];
            
            if (parms["pin"] > 0) {
                for (int i = 0; i < parms["pin"]; ++i)
                {
                    terms.clear();
                    terms.push_back( std::make_pair(i, h0*c[(i+2*(i/W))%3]) );
                    mm.add_term(terms);
                }
            }
            
            if (parms["pin"] < 0) {
                for (int i = 0; i < -parms["pin"]; ++i)
                {
                    terms.clear();
                    int t = (i % 2 == 0 ? 0 : 2);
                    terms.push_back( std::make_pair(i, h0*c[t]) );
                    mm.add_term(terms);
                }
            }
        }
            
    private:
        double Jbl, Jbq;
        op_pairs ops;
        BaseParameters & parms;
    };
}

#endif
