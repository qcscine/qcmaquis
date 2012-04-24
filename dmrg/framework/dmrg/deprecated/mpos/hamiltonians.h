/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef HAMILTONIANS_H
#define HAMILTONIANS_H

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/deprecated/mpos/generate_mpo.h"

#include <vector>
#include <utility>

namespace mpos {
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
    class Spin1BlBq : public Hamiltonian<Matrix, U1>
    {
    public:
        typedef typename Hamiltonian<Matrix, U1>::op_t op_t;
        typedef std::pair<op_t, op_t> op_pair;
        typedef std::vector<op_pair> op_pairs;
        
        Spin1BlBq(double Jbl_, double Jbq_, double h0_) : Jbl(Jbl_), Jbq(Jbq_), h0(h0_)
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
        
        void push_extra_terms(MPOMaker<Matrix, U1> & mm)
        {
            maquis::cout << "Adding extra term." << std::endl;
            std::vector<std::pair<std::size_t, op_t> > terms;
            block_matrix<Matrix, U1> sz, szz;
            
            sz.insert_block(Matrix(1, 1, 1), 1, 1);
            sz.insert_block(Matrix(1, 1, 0), 0, 0);
            sz.insert_block(Matrix(1, 1, -1), -1, -1);
            
            gemm(sz, sz, szz);
            
            terms.push_back( std::make_pair(1, h0*sz) );
            
            mm.add_term(terms);
        }
            
    private:
        double Jbl, Jbq, h0;
        op_pairs ops;
    };
}

#endif
