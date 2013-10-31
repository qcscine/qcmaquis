/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MPO_OPS_H
#define MPO_OPS_H

#include "mpo.h"
#include "mpotensor.h"

template<class Matrix, class SymmGroup>
std::string identify_op(block_matrix<Matrix, SymmGroup> const & m)
{
    if (m.n_blocks() == 2)
        return "I";
    else {
        typename SymmGroup::charge c1 = m.left_basis()[0].first, c2 = m.right_basis()[0].first;
        if (c1 == 1 && c2 == 0)
            return "c";
        else if (c1 == 0 && c2 == 1)
            return "d";
        else if (c1 == 1 && c2 == 1)
            return "cd";
        else if (c1 == 0 && c2 == 0)
            return "dc";
    }
    
    return "wtf?";
}

template<class Matrix, class SymmGroup>
void follow_mpo(MPO<Matrix, SymmGroup> const & mpo,
                std::string s = std::string(),
                int p = 0, int start = 0)
{
    for (size_t k = 0; k < mpo[p].col_dim(); ++k)
    {
        if (mpo[p].at(start,k).op.n_blocks() == 0)
            continue;
        
        std::ostringstream oss;
//        oss << mpo[p](start, k) << std::endl;
//        oss << "(" << start << "," << k << ") ";
        oss << " " << identify_op(mpo[p].at(start, k).op) << " ";
        if (p+1 < mpo.length())
            follow_mpo(mpo, s+oss.str(), p+1, k);
        else
            maquis::cout << s+oss.str() << std::endl;
    }
}

template<class Matrix, class SymmGroup>
void follow_and_print_terms(MPO<Matrix, SymmGroup> const& mpo, int p, int b1, int b2, std::string s="", typename MPOTensor<Matrix,SymmGroup>::value_type scale=1.)
{
    std::stringstream ss;
    ss << s;
    
    if (p > -1) {
        MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup> access = mpo[p].at(b1,b2);
        scale *= access.scale;
        ss << " {" << mpo[p].tag_number(b1,b2) << "}(" << p << ")";
    }
    
    if (p == mpo.size()-1) {
        maquis::cout << "---" << std::endl;
        maquis::cout << "scale: " << scale << std::endl;
        maquis::cout << "term: "  << ss.str() << std::endl;
        return;
    }
    
    typedef typename MPOTensor<Matrix, SymmGroup>::row_proxy row_proxy;
    typedef typename MPOTensor<Matrix, SymmGroup>::col_proxy col_proxy;
    row_proxy myrow = mpo[p+1].row(b2);
    for (typename row_proxy::const_iterator row_it = myrow.begin(); row_it != myrow.end(); ++row_it)
        follow_and_print_terms(mpo, p+1, b2, row_it.index(), ss.str(), scale);
}


template<class Matrix, class SymmGroup>
void cleanup_mpo_(MPO<Matrix, SymmGroup> const & in_mpo,
                  MPO<Matrix, SymmGroup> & out_mpo,
                  std::vector<boost::tuple<int, int, block_matrix<Matrix, SymmGroup> > > & ops,
                  int p, int start)
{
    for (std::size_t k = 0; k < in_mpo[p].col_dim(); ++k)
    {
        if (!in_mpo[p].has(start,k))
            continue;
        if (in_mpo[p].at(start,k).op.n_blocks() == 0)
            continue;
        
        ops[p] = boost::make_tuple(start, k, in_mpo[p].at(start, k).op * in_mpo[p].at(start, k).scale);
        
        if (p+1 < in_mpo.length())
            cleanup_mpo_(in_mpo, out_mpo, ops, p+1, k);
        else
        {
            assert( ops.size() == out_mpo.length() );
            for (std::size_t t = 0; t < in_mpo.length(); ++t) {
                //block_matrix<Matrix, SymmGroup> & out_b = out_mpo[t](boost::tuples::get<0>(ops[t]),
                //                                                     boost::tuples::get<1>(ops[t]));
                
                //if (out_b.n_blocks() == 0)
                //    out_b = boost::tuples::get<2>(ops[t]);
                if (out_mpo[t].at(boost::tuples::get<0>(ops[t]), boost::tuples::get<1>(ops[t])).op.n_blocks() == 0)
                    out_mpo[t].set(boost::tuples::get<0>(ops[t]), boost::tuples::get<1>(ops[t]), boost::tuples::get<2>(ops[t]));
            }
        }
    }
}

template<class Matrix, class SymmGroup>
MPO<Matrix, SymmGroup> cleanup_mpo(MPO<Matrix, SymmGroup> const & mpo)
{
    MPO<Matrix, SymmGroup> ret(mpo.length());
    for (std::size_t p = 0; p < ret.length(); ++p)
        ret[p] = MPOTensor<Matrix, SymmGroup>(mpo[p].row_dim(), mpo[p].col_dim());
    
    std::vector<boost::tuple<int, int, block_matrix<Matrix, SymmGroup> > > prempo(mpo.length());
    cleanup_mpo_(mpo, ret, prempo, 0, 0);
    return ret;
}

template<class Matrix, class SymmGroup>
MPO<Matrix, SymmGroup>
square_mpo(MPO<Matrix, SymmGroup> const & mpo)
{
    typedef typename SymmGroup::charge charge;
    typedef typename MPOTensor<Matrix, SymmGroup>::row_proxy row_proxy;
    typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
    
    size_t L = mpo.length();
    
    MPO<Matrix, SymmGroup> sq(L);
    
    for (size_t p = 0; p < L; ++p)
    {
        MPOTensor<Matrix, SymmGroup> const & inp = mpo[p];
        maquis::cout << "MPOTensor " << inp.row_dim()*inp.row_dim() << " " << inp.col_dim()*inp.col_dim() << std::endl;
        MPOTensor<Matrix, SymmGroup> ret(inp.row_dim()*inp.row_dim(),
                                         inp.col_dim()*inp.col_dim());
        
        for (index_type r1 = 0; r1 < inp.row_dim(); ++r1)
        {
            row_proxy row1 = inp.row(r1);
            for (index_type r2 = 0; r2 < inp.row_dim(); ++r2)
            {
                row_proxy row2 = inp.row(r2);
                for (typename row_proxy::const_iterator it1 = row1.begin(); it1 != row1.end(); ++it1)
                {
                    index_type c1 = it1.index();
                    for (typename row_proxy::const_iterator it2 = row2.begin(); it2 != row2.end(); ++it2) {
                        index_type c2 = it2.index();

                        assert(inp.has(r1, c1));
                        assert(inp.has(r2, c2));
                        
                        block_matrix<Matrix, SymmGroup> t;
                        gemm(inp.at(r1, c1).op, inp.at(r2, c2).op, t);
                        if (t.n_blocks() > 0)
                            ret.set(r1*inp.row_dim()+r2, c1*inp.col_dim()+c2, 
                                        t * (inp.at(r1, c1).scale * inp.at(r2, c2).scale));
                    }
                }
            }
        }
        
        sq[p] = ret;
    }
    
    maquis::cout << "Done squaring." << std::endl;
    
    sq = cleanup_mpo(sq);
    maquis::cout << "Done cleaning up." << std::endl;

    return sq;
}

template<class Matrix, class SymmGroup>
MPO<Matrix, SymmGroup>
zero_after(MPO<Matrix, SymmGroup> mpo, int p0)
{
    typedef typename MPOTensor<Matrix, SymmGroup>::CSRMatrix CSRMatrix;
    typedef typename MPOTensor<Matrix, SymmGroup>::CSCMatrix CSCMatrix;

    maquis::cout << "Zeroing out MPO after site " << p0 << std::endl;

    for (int p = p0+1; p < mpo.size(); ++p) {
        for (int k = 2; k < mpo[p].row_dim(); ++k)
            for (int l = 2; l < mpo[p].col_dim(); ++l)
                if (mpo[p].has(k,l))
                    mpo[p].set(k,l, mpo[p].at(k,l).op, 0.0);
    
        if (mpo[p].has(0,1))
            mpo[p].set(0,1, mpo[p].at(0,1).op, 0.0);
    }
    
    return mpo;
}

#endif
