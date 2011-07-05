/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/
 
#ifndef MPO_H
#define MPO_H

#include <vector>
#include <set>

#include "mp_tensors/mpotensor.h"

template<class Matrix, class SymmGroup>
class MPO : public std::vector<MPOTensor<Matrix, SymmGroup> >
{
public:
    typedef MPOTensor<Matrix, SymmGroup> elem_type;
    
    MPO(std::size_t L, elem_type elem = elem_type())
    : std::vector<elem_type>(L, elem)
    { }
    
    std::size_t length() const { return this->size(); }
    
    void compress(double cutoff)
    {
        calc_charges();
        
        for (int p = 0; p < this->size()-1; ++p) {
            MPOTensor<Matrix, SymmGroup> b1 = (*this)[p], b2 = (*this)[p+1];
            
            block_matrix<Matrix, SymmGroup> left = make_left_matrix(p);
            block_matrix<Matrix, SymmGroup> right = make_right_matrix(p+1);
            
            block_matrix<Matrix, SymmGroup> M, U, V;

            block_matrix<typename blas::associated_diagonal_matrix<Matrix>::type, SymmGroup> S, Sqrt;
            gemm(left, right, M);
            svd_truncate(M, U, V, S, cutoff, 100000, false);
            Sqrt = sqrt(S);
            gemm(U, Sqrt, left);
            gemm(Sqrt, V, right);
            
            zout << "MPO bond truncation: " << bond_indices[p+1].sum_of_sizes() << " -> ";
            replace_pair(left, right, p);
            zout << bond_indices[p+1].sum_of_sizes() << std::endl;
        }
    }
    
private:
    std::vector<std::map<std::size_t, typename SymmGroup::charge> > bond_index_charges;
    std::vector<Index<SymmGroup> > bond_indices;
    
    void calc_charges()
    {
        using std::size_t;
        
        size_t L = this->size();
        bond_index_charges.clear();
        bond_index_charges.resize(L+1);
        
        bond_index_charges[0][0] = SymmGroup::SingletCharge;
        
        for (size_t p = 1; p <= L; ++p)
        {
            for (size_t c = 0; c < (*this)[p-1].col_dim(); ++c) {
                std::set<typename SymmGroup::charge> charge_diffs;
                for (size_t r = 0; r < (*this)[p-1].row_dim(); ++r) {
                    assert( bond_index_charges[p-1].count(r) > 0 );
                    if (!(*this)[p-1].has(r,c))
                        continue;
                    for (size_t b = 0; b < (*this)[p-1](r, c).n_blocks(); ++b) {
                        charge_diffs.insert(SymmGroup::fuse(bond_index_charges[p-1][r],
                                                            SymmGroup::fuse((*this)[p-1](r,c).left_basis()[b].first,
                                                                            -(*this)[p-1](r,c).right_basis()[b].first)));
//                        std::cout << r << " " << c << std::endl;
//                        std::cout << bond_index_charges[p-1][r] << std::endl;
//                        std::cout << (*this)[p-1](r,c).left_basis()[b].first << std::endl;
//                        std::cout << (*this)[p-1](r,c).right_basis()[b].first << std::endl;
//                        std::copy(charge_diffs.begin(), charge_diffs.end(),
//                                  std::ostream_iterator<typename SymmGroup::charge>(std::cout, " "));
//                        std::cout << std::endl << std::endl;
                    }
                }
#ifndef NDEBUG
                if (charge_diffs.size() > 1) {
                    std::copy(charge_diffs.begin(), charge_diffs.end(),
                              std::ostream_iterator<typename SymmGroup::charge>(std::cout, " "));
                    zout << std::endl;
                }
                assert( charge_diffs.size() <= 1 );
#endif
                if (charge_diffs.size() == 1)
                    bond_index_charges[p][c] = *charge_diffs.begin();
                else
                    bond_index_charges[p][c] = SymmGroup::SingletCharge; //bond_index_charges[p-1][c];
            }
        }
        
        bond_indices.clear();
        bond_indices.resize(L+1);
        for (size_t p = 0; p <= L; ++p)
        {
            Index<SymmGroup> & index = bond_indices[p];
            for (typename std::map<std::size_t, typename SymmGroup::charge>::iterator it
                 = bond_index_charges[p].begin();
                 it != bond_index_charges[p].end();
                 ++it)
                if (index.has(it->second))
                    index[index.position(it->second)] = std::make_pair(it->second, index.size_of_block(it->second)+1);
                else
                    index.insert(std::make_pair(it->second, 1));
            zout << index << std::endl;
        }
    }
    
    block_matrix<Matrix, SymmGroup> make_left_matrix(std::size_t p)
    {
        using std::size_t;
        typedef typename SymmGroup::charge charge;
        
        Index<SymmGroup> phys_i;
        for (size_t r = 0; r < (*this)[p].row_dim(); ++r)
            for (size_t c = 0; c < (*this)[p].col_dim(); ++c)
            {
                if (!(*this)[p].has(r,c))
                    continue;
                for (size_t cs = 0; cs < (*this)[p](r, c).left_basis().size(); ++cs) {
                    std::pair<charge, size_t> sector = (*this)[p](r, c).left_basis()[cs];
                    if (! phys_i.has(sector.first))
                        phys_i.insert(sector);
                }
            }
        
        Index<SymmGroup> left_i = phys_i * adjoin(phys_i) * bond_indices[p];
        Index<SymmGroup> right_i = bond_indices[p+1];
        
        left_i = common_subset(left_i, right_i);
        
        block_matrix<Matrix, SymmGroup> ret(left_i, right_i);
        
        std::map<charge, size_t> visited_c_basis;
        for (size_t c = 0; c < (*this)[p].col_dim(); ++c) {
            int outr = -1;
            for (size_t r = 0; r < (*this)[p].row_dim(); ++r)
                for (size_t ls = 0; ls < phys_i.size(); ++ls)
                    for (size_t rs = 0; rs < phys_i.size(); ++rs) {
                        charge lc = SymmGroup::fuse(bond_index_charges[p][r],
                                                    SymmGroup::fuse(phys_i[ls].first,
                                                                    -phys_i[rs].first));
                        charge rc = bond_index_charges[p+1][c];
                        
                        if (lc != rc)
                            continue;
                        
                        outr++;
                        
                        if (! (*this)[p].has(r,c))
                            continue;
                        if (! (*this)[p](r,c).has_block(phys_i[ls].first, phys_i[rs].first) )
                            continue;                       
                        
                        std::size_t cs = (*this)[p](r, c).left_basis().position(phys_i[ls].first);
                        
                        assert( lc == rc );
                        assert( outr < left_i.size_of_block(lc) );
                        
                        // ...for now...
                        assert( num_rows((*this)[p](r,c)[cs]) == 1 );
                        assert( num_cols((*this)[p](r,c)[cs]) == 1 );
                        
                        //ret(std::make_pair(lc, outr),
                        //    std::make_pair(rc, visited_c_basis[rc]));
                        //(*this)[p](r,c)[cs](0,0);
                        // ambient_write_only //
                        ret(std::make_pair(lc, outr),
                            std::make_pair(rc, visited_c_basis[rc])) =
                        (*this)[p](r,c)[cs](0,0);
                        
//                        std::cout << (*this)[p](r,c)[cs](0,0) << " | ";
//                        std::cout << r << " " << c << " " << phys_i[ls].first << " " << phys_i[rs].first;
//                        std::cout << " -> ";
//                        std::cout << "(" << lc << "," << outr << ") (" << rc << "," << visited_c_basis[rc] << ")" << std::endl;
                    }
            visited_c_basis[bond_index_charges[p+1][c]]++;
        }
        
//        std::cout << ret << std::endl;
//        std::cout << "###" << std::endl;
        
        return ret;
    }
    
    block_matrix<Matrix, SymmGroup> make_right_matrix(std::size_t p)
    {
        using std::size_t;
        typedef typename SymmGroup::charge charge;
        
        Index<SymmGroup> phys_i;
        for (size_t r = 0; r < (*this)[p].row_dim(); ++r)
            for (size_t c = 0; c < (*this)[p].col_dim(); ++c)
            {
                if (!(*this)[p].has(r,c))
                    continue;
                for (size_t cs = 0; cs < (*this)[p](r, c).left_basis().size(); ++cs) {
                    std::pair<charge, size_t> sector = (*this)[p](r, c).left_basis()[cs];
                    if (! phys_i.has(sector.first))
                        phys_i.insert(sector);
                }
            }
        
        Index<SymmGroup> left_i = bond_indices[p];
        Index<SymmGroup> right_i = adjoin(phys_i) * phys_i * bond_indices[p+1];
        
        right_i = common_subset(right_i, left_i);
        
        block_matrix<Matrix, SymmGroup> ret(left_i, right_i);
        
        std::map<charge, size_t> visited_r_basis;
        for (size_t r = 0; r < (*this)[p].row_dim(); ++r) {
            int outc = -1;
            for (size_t c = 0; c < (*this)[p].col_dim(); ++c)
                for (size_t ls = 0; ls < phys_i.size(); ++ls)
                    for (size_t rs = 0; rs < phys_i.size(); ++rs) {
                        charge lc = bond_index_charges[p][r];
                        charge rc = SymmGroup::fuse(bond_index_charges[p+1][c],
                                                    SymmGroup::fuse(-phys_i[ls].first,
                                                                    phys_i[rs].first));
                        
                        if (lc != rc)
                            continue;
                        
                        outc++;
                        
                        if (! (*this)[p].has(r,c))
                            continue;
                        if (! (*this)[p](r, c).has_block(phys_i[ls].first, phys_i[rs].first) )
                            continue;
                        
                        std::size_t cs = (*this)[p](r, c).left_basis().position(phys_i[ls].first);
                        
                        assert( lc == rc );
                        assert( outc < right_i.size_of_block(rc) );
                        
                        // ...for now...
                        assert( num_rows((*this)[p](r,c)[cs]) == 1 );
                        assert( num_cols((*this)[p](r,c)[cs]) == 1 );
                        
                        // ambient_write_only //
                        ret(std::make_pair(lc, visited_r_basis[lc]),
                            std::make_pair(rc, outc)) =
                        (*this)[p](r,c)[cs](0,0);
                        
//                        cout << (*this)[p](r,c)[cs](0,0) << " | ";
//                        cout << r << " " << c << " " << phys_i[ls].first << " " << phys_i[rs].first;
//                        cout << " -> ";
//                        cout << "(" << lc << "," << visited_r_basis[lc] << ") (" << rc << "," << outc << ")" << endl;
                    }
            visited_r_basis[bond_index_charges[p][r]]++;
        }
        
//        cout << ret << endl;
//        cout << "###" << endl;
        
        return ret;
    }
    
    void replace_pair(block_matrix<Matrix, SymmGroup> & left,
                      block_matrix<Matrix, SymmGroup> & right,
                      std::size_t p)
    {
        using std::size_t;
        typedef typename SymmGroup::charge charge;
        
        Index<SymmGroup> phys_i;
        for (size_t r = 0; r < (*this)[p].row_dim(); ++r)
            for (size_t c = 0; c < (*this)[p].col_dim(); ++c)
            {
                for (size_t cs = 0; cs < (*this)[p](r, c).left_basis().size(); ++cs) {
                    std::pair<charge, size_t> sector = (*this)[p](r, c).left_basis()[cs];
                    if (! phys_i.has(sector.first))
                        phys_i.insert(sector);
                }
            }
        
        assert( left.right_basis() == right.left_basis() );
        bond_indices[p+1] = left.right_basis();
        {
            std::size_t count = 0;
            bond_index_charges[p+1].clear();
            for (typename Index<SymmGroup>::basis_iterator it = left.right_basis().basis_begin();
                 !it.end(); ++it)
                bond_index_charges[p+1][count++] = (*it).first;
        }
        
        (*this)[p] = MPOTensor<Matrix, SymmGroup>((*this)[p].row_dim(),
                                                  bond_indices[p+1].sum_of_sizes());
        
        std::map<charge, size_t> visited_c_basis;
        for (size_t c = 0; c < (*this)[p].col_dim(); ++c) {
            int outr = -1;
            for (size_t r = 0; r < (*this)[p].row_dim(); ++r)
                for (size_t ls = 0; ls < phys_i.size(); ++ls)
                    for (size_t rs = 0; rs < phys_i.size(); ++rs)
                    {
                        charge lc = SymmGroup::fuse(bond_index_charges[p][r],
                                                    SymmGroup::fuse(phys_i[ls].first,
                                                                    -phys_i[rs].first));
                        charge rc = bond_index_charges[p+1][c];
                        
                        if (lc != rc)
                            continue;
                        
                        outr++;
                        
                        typename Matrix::value_type val = left(std::make_pair(lc, outr),
                                                               std::make_pair(rc, visited_c_basis[rc]));
                        
                        if (fabs(val) > 1e-40) {
                            block_matrix<Matrix, SymmGroup> & block = (*this)[p](r,c);
                            charge blc = phys_i[ls].first, brc = phys_i[rs].first;
                            
                            block.insert_block(Matrix(1, 1, val), blc, brc);
                            
//                            cout << val << " | ";
//                            cout << r << " " << c << " " << phys_i[ls].first << " " << phys_i[rs].first;
//                            cout << " <- ";
//                            cout << "(" << lc << "," << outr << ") (" << rc << "," << visited_c_basis[rc] << ")" << endl;
                        }
                    }
            visited_c_basis[bond_index_charges[p+1][c]]++;
        }
        
        (*this)[p+1] = MPOTensor<Matrix, SymmGroup>(bond_indices[p+1].sum_of_sizes(),
                                                    (*this)[p+1].col_dim());
        
        std::map<charge, size_t> visited_r_basis;
        for (size_t r = 0; r < (*this)[p+1].row_dim(); ++r) {
            int outc = -1;
            for (size_t c = 0; c < (*this)[p+1].col_dim(); ++c)
                for (size_t ls = 0; ls < phys_i.size(); ++ls)
                    for (size_t rs = 0; rs < phys_i.size(); ++rs)
                    {
                        charge lc = bond_index_charges[p+1][r];
                        charge rc = SymmGroup::fuse(bond_index_charges[p+2][c],
                                                    SymmGroup::fuse(-phys_i[ls].first,
                                                                    phys_i[rs].first));
                        
                        if (lc != rc)
                            continue;
                        
                        outc++;
                        
                        typename Matrix::value_type val = right(std::make_pair(lc, visited_r_basis[lc]),
                                                                std::make_pair(rc, outc));
                        
                        if (fabs(val) > 1e-40) {
                            block_matrix<Matrix, SymmGroup> & block = (*this)[p+1](r,c);
                            charge blc = phys_i[ls].first, brc = phys_i[rs].first;
                            
                            block.insert_block(Matrix(1, 1, val), blc, brc);
                            
//                            cout << val << " | ";
//                            cout << r << " " << c << " " << phys_i[ls].first << " " << phys_i[rs].first;
//                            cout << " <- ";
//                            cout << "(" << lc << "," << visited_r_basis[lc] << ") (" << rc << "," << outc << ")" << endl;
                        }
                    }
            visited_r_basis[bond_index_charges[p+1][r]]++;
        }
    }
};

#endif
