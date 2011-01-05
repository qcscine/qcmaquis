#ifndef HAMILTONIANS_H
#define HAMILTONIANS_H

#include "block_matrix.h"

#include <boost/tuple/tuple.hpp>

#include <vector>
#include <utility>

template<class Matrix, class SymmGroup>
class Hamiltonian
{
public:
    typedef block_matrix<Matrix, SymmGroup> op_t;
    typedef boost::tuple<size_t, size_t, op_t> block;
    typedef std::vector<
	    std::pair<
    		block_matrix<Matrix, SymmGroup>,
		    block_matrix<Matrix, SymmGroup>
	    >
    > op_pairs;
    
    virtual op_pairs get_ops() = 0;
};

template<class Matrix>
class Heisenberg : public Hamiltonian<Matrix, U1>
{
public:
    Heisenberg(double Jxy_, double Jz_) : Jxy(Jxy_), Jz(Jz_) { }
    
    typedef typename Hamiltonian<Matrix, U1>::op_pairs op_pairs;
    
    op_pairs get_ops()
    {
        using boost::tuples::make_tuple;
        
        block_matrix<Matrix, U1> ident, splus, sminus, sz;
        
        ident.insert_block( make_tuple(Matrix(1, 1, 1), -1, -1) );
        ident.insert_block( make_tuple(Matrix(1, 1, 1), 1, 1) );
        
        splus.insert_block( make_tuple(Matrix(1, 1, 1), -1, 1) );
        
        sminus.insert_block( make_tuple(Matrix(1, 1, 1), 1, -1) );
        
        sz.insert_block( make_tuple(Matrix(1, 1, 0.5), 1, 1) );
        sz.insert_block( make_tuple(Matrix(1, 1, -0.5), -1, -1) );
        
        op_pairs ret;
        
        ret.push_back(make_pair(ident, ident));
        ret.push_back(make_pair(Jxy/2*splus, sminus));
        ret.push_back(make_pair(Jxy/2*sminus, splus));
        ret.push_back(make_pair(Jz*sz, sz));
        
        return ret;
    }
    
private:
    double Jxy, Jz;
};

template<class Matrix>
class Spin1BlBq : public Hamiltonian<Matrix, U1>
{
public:
    Spin1BlBq(double Jbl_, double Jbq_) : Jbl(Jbl_), Jbq(Jbq_) { }
    
    typedef typename Hamiltonian<Matrix, U1>::op_pairs op_pairs;
    
    op_pairs get_ops()
    {
        using boost::tuples::make_tuple;
        
        block_matrix<Matrix, U1> ident, splus, sminus, sz, spp, smm, spm, smp, szz, szp, spz, szm, smz;
        
        ident.insert_block( make_tuple(Matrix(1, 1, 1), -1, -1) );
        ident.insert_block( make_tuple(Matrix(1, 1, 1), 0, 0) );
        ident.insert_block( make_tuple(Matrix(1, 1, 1), 1, 1) );
        
        splus.insert_block( make_tuple(Matrix(1, 1, 1), -1, 0) );
        splus.insert_block( make_tuple(Matrix(1, 1, 1), 0, 1) );
        
        sminus.insert_block( make_tuple(Matrix(1, 1, 1), 1, 0) );
        sminus.insert_block( make_tuple(Matrix(1, 1, 1), 0, -1) );
        
        sz.insert_block( make_tuple(Matrix(1, 1, 1), 1, 1) );
        sz.insert_block( make_tuple(Matrix(1, 1, 0), 0, 0) );
        sz.insert_block( make_tuple(Matrix(1, 1, -1), -1, -1) );
        
        gemm(splus, splus, spp);
        gemm(sminus, sminus, smm);
        gemm(splus, sminus, spm);
        gemm(sminus, splus, smp);
        gemm(sz, sz, szz);
        gemm(sz, splus, szp);
        gemm(splus, sz, spz);
        gemm(sz, sminus, szm);
        gemm(sminus, sz, smz);
        
        op_pairs ret;
        
#define term(a,b) ret.push_back(make_pair(a, b))
        
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
        
        return ret;
    }

private:
    double Jbl, Jbq;
};

#endif
