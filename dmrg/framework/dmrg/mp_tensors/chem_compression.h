/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef FAST_COMPRESSION_H
#define FAST_COMPRESSION_H

#include <unistd.h>
#include <vector>
#include <functional>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"

#include "dmrg/models/lattice.h"

#include "dmrg/mp_tensors/mpotensor.h"


template <class SymmGroup>
double mem_index(Index<SymmGroup> const & li, Index<SymmGroup> const & ri)
{
    typedef typename SymmGroup::charge charge;
    typedef std::size_t size_type;
    std::vector<size_type> lefts = li.sizes(), rights = ri.sizes();
    size_type foot_print = 0;
    for(int s=0; s < lefts.size(); s++) {
        foot_print += lefts[s]*rights[s];
        maquis::cout //<< ((lefts[s] > 10000) ? (boost::lexical_cast<std::string>(lefts[s]/1000) + std::string("K"))
                     //   : boost::lexical_cast<std::string>(lefts[s]))
                     << lefts[s]
                     << "x"
                     //<< ((rights[s] > 10000) ? (boost::lexical_cast<std::string>(rights[s]/1000) + std::string("K"))
                     //   : boost::lexical_cast<std::string>(rights[s]))
                     << rights[s]
                     << "("
                     << ((lefts[s]*rights[s] > 10000) ? (boost::lexical_cast<std::string>(lefts[s]*rights[s]/1000) + std::string("K"))
                        : boost::lexical_cast<std::string>(lefts[s]*rights[s]))
                     << ") ";
    }

    double fp = double(foot_print)*8./double(1u<<30);
    maquis::cout << "footprint: " << fp << " GB" << std::endl;
    return fp;
}

template <class Matrix, class SymmGroup>
void mem_footprint(block_matrix<Matrix, SymmGroup> const & rval)
{
    mem_index<SymmGroup>(rval.left_basis(), rval.right_basis());
}

template <class SymmGroup>
void mem_footprint(std::vector<Index<SymmGroup> > const & bi, Index<SymmGroup> const & phys)
{
    for(int p=0; p < bi.size()-2; ++p) {
        Index<SymmGroup> left_i = phys*adjoin(phys)*bi[p];
        Index<SymmGroup> mid_i  = bi[p+1];
        left_i = common_subset(left_i, mid_i);
        Index<SymmGroup> right_i  = adjoin(phys)*phys*bi[p+2];
        right_i = common_subset(mid_i, right_i);

        maquis::cout << "L" << p << " ";
        mem_index(left_i, mid_i);
        maquis::cout << "R" << p+1 << " ";
        mem_index(mid_i, right_i);
    }
}

namespace compressor_detail
{
    template <class SymmGroup>
    struct charge_pair_gt
    {
        typedef typename SymmGroup::charge charge;
        typedef unsigned charge_cnt_t;

        bool operator()(std::pair<charge, charge> const & i,
                        std::pair<charge, charge> const & j) const
        {
            using namespace boost::tuples;
            if (i.first > j.first)
                return true;
            else if (i.first < j.first)
                return false;
            else
                return i.second > j.second;
        }
    };

    struct pair_cmp_inv
    {  
        bool operator()(std::pair<std::size_t, std::size_t> const & i,
                        std::pair<std::size_t, std::size_t> const & j) const
        {  
            if (i.second < j.second)
                return true;
            else if (i.second > j.second)
                return false;
            else
                return i.first < j.first;
        }
    };
}

template<class Matrix, class SymmGroup>
class charge_overlap
{
    typedef typename SymmGroup::charge charge;

public:
    charge_overlap() {}
    charge_overlap(Index<SymmGroup> const & phys)
    {
        throw std::runtime_error("Charge overlap only implemented for TwoU1 symmetry");
    }

    std::size_t operator()(charge c1, charge c2) { return 0; }
    std::vector<charge> get_charge_deltas() const { return std::vector<charge>(); }

};

template<class Matrix>
class charge_overlap<Matrix, TwoU1>
{
    typedef typename TwoU1::charge charge;
    typedef typename compressor_detail::charge_pair_gt<TwoU1>::charge_cnt_t charge_cnt_t;
    typedef typename Index<TwoU1>::value_type index_value_t;
    typedef typename Index<TwoU1>::iterator index_iterator;

public:
    charge_overlap(Index<TwoU1> const & phys) : phys_i(phys)
    {
        Index<TwoU1> one_site_mutations = phys_i * adjoin(phys_i);
        Index<TwoU1> two_site_mutations = one_site_mutations * phys_i * adjoin(phys_i);
        // TwoU1 only
        index_iterator new_end =
            std::remove_if(two_site_mutations.begin(), two_site_mutations.end(), boost::bind(&index_value_t::second, boost::lambda::_1) < 6);
        two_site_mutations.erase(new_end, two_site_mutations.end());

        maquis::cout << "One site and twosite mutations:\n";
        maquis::cout << one_site_mutations << std::endl;
        maquis::cout << two_site_mutations << std::endl;

        typedef typename std::vector<charge>::iterator charge_vec_it;
        std::vector<charge> ss_mut_vec = one_site_mutations.charges();
        std::vector<charge> ts_mut_vec = two_site_mutations.charges();
        charge_deltas = ts_mut_vec;
        charge_vec_it it1, it2, ss_it;

        // build (global) lookup table
        for( it1 = ts_mut_vec.begin(); it1 != ts_mut_vec.end(); ++it1)
            for( it2 = ts_mut_vec.begin(); it2 != ts_mut_vec.end(); ++it2)
            {
                charge_cnt_t outr = 0;

                for (size_t ls = 0; ls < phys_i.size(); ++ls)
                    for (size_t rs = 0; rs < phys_i.size(); ++rs) {
                        charge lc = TwoU1::fuse(*it1, TwoU1::fuse(phys_i[ls].first, -phys_i[rs].first));
                        if (lc == *it2) outr++;
                    }

                op_chain_map_hash[std::make_pair(*it1, *it2)] = outr;
                op_chain_map[std::make_pair(*it1, *it2)] = outr;
            }

        // check & print lookup table
        maquis::cout << std::endl;
        for( it1 = ts_mut_vec.begin(); it1 != ts_mut_vec.end(); ++it1) {
            for( it2 = ts_mut_vec.begin(); it2 != ts_mut_vec.end(); ++it2) {
                if( op_chain_map_hash[std::make_pair(*it1, *it2)] != op_chain_map_hash[std::make_pair(*it2, *it1)] )
                    exit(1);
                maquis::cout << op_chain_map_hash[std::make_pair(*it1, *it2)] << " ";
            }
            maquis::cout << std::endl;
        }

    }

    std::size_t operator()(charge c1, charge c2)
    {
        return op_chain_map_hash[std::make_pair(c1,c2)];
        //return op_chain_map[std::make_pair(c1,c2)];
    }

    std::vector<charge> get_charge_deltas() const { return charge_deltas; }
    
private:
    Index<TwoU1> phys_i;
    std::map< std::pair<charge, charge>, charge_cnt_t, compressor_detail::charge_pair_gt<TwoU1> > op_chain_map;
    boost::unordered_map<std::pair<charge, charge>, charge_cnt_t> op_chain_map_hash;
    std::vector<charge> charge_deltas;
};

template <class Matrix, class SymmGroup>
class MPOIndexer : public std::vector< std::vector<typename SymmGroup::charge> >
{
public:
    typedef std::vector<typename SymmGroup::charge> element_type;
    typedef typename MPO<Matrix, SymmGroup>::elem_type::data_t::const_iterator map_it_t;
    typedef unsigned mpo_index_t;
    typedef int offset_t;

    MPOIndexer() {}
    
    MPOIndexer(MPO<Matrix, SymmGroup> const & mpo_in)
    {
        mpo_index_t L = mpo_in.size();
        this->resize(L+1);
        bond_indices.resize(L+1);

        (*this)[0].resize(1);
        (*this)[L].resize(1);
        for ( int k=1; k < (*this).size()-1; ++k)
            (*this)[k].resize(mpo_in[k-1].col_dim());
        
        (*this)[0][0] = SymmGroup::IdentityCharge;

        Timer ctim("MPOIndexer ctor");
        ctim.begin();
        
        for (mpo_index_t p = 1; p <= L; ++p)
        {
            // loop over all operator terms that start here(p-1) (row_dim == 0)
            for (map_it_t it=mpo_in[p-1].data_.lower_bound(std::make_pair(0,0));
                 it != mpo_in[p-1].data_.lower_bound(std::make_pair(1,0)); ++it){

                mpo_index_t r = 0, c = 0;
                offset_t tp = p-2;

                do {
                    // go to next term
                    ++tp;
                    r = c;
                    if (r > 0) {
                        // there should be only one element until lower_bound(pair(r+1,))
                        map_it_t next_term = mpo_in[tp].data_.lower_bound(std::make_pair(r, 0));
                        c = next_term->first.second;
                        assert( next_term->first.first == r);
                    } else c = it->first.second;

                    typename SymmGroup::charge charge_diff = 
                                SymmGroup::fuse((*this)[tp][r],
                                SymmGroup::fuse(mpo_in[tp](r,c).left_basis()[0].first,
                                                -mpo_in[tp](r,c).right_basis()[0].first));

                    (*this)[tp+1][c] = charge_diff;
                // while the operator chain starting here(p-1) goes on, follow the terms
                } while (c > 1);
            }
        }
        
        bond_indices.clear();
        bond_indices.resize(L+1);
        for (mpo_index_t p = 0; p <= L; ++p)
        {
            Index<SymmGroup> & index = bond_indices[p];
            for (typename element_type::iterator it= (*this)[p].begin(); it != (*this)[p].end(); ++it)
                if (index.has(*it))
                    index[index.position(*it)] = std::make_pair(*it, index.size_of_block(*it)+1);
                else
                    index.insert(std::make_pair(*it, 1));
        }
        ctim.end();
    }

    void update_bond(Index<SymmGroup> const & rhs, mpo_index_t p)
    {
        bond_indices[p] = rhs;
        mpo_index_t count = 0;
        (*this)[p].clear();
        for (typename Index<SymmGroup>::basis_iterator it = rhs.basis_begin(); !it.end(); ++it)
            (*this)[p][count++] = (*it).first;
    }

    Index<SymmGroup> & index(mpo_index_t p) { return bond_indices[p]; } 
    Index<SymmGroup> const & index(mpo_index_t p) const { return bond_indices[p]; } 

private:
    std::vector< Index<SymmGroup> > bond_indices;
};

///////////////////////////////////////////////////////////////////////////////

template<class Matrix, class SymmGroup>
class column_iterator
{
    typedef typename SymmGroup::charge charge;
    typedef typename MPOIndexer<Matrix, SymmGroup>::element_type bond_charge_t;
    typedef typename MPOIndexer<Matrix, SymmGroup>::mpo_index_t mpo_index_t;
    typedef typename MPOIndexer<Matrix, SymmGroup>::offset_t offset_t;
    typedef typename MPOIndexer<Matrix, SymmGroup>::map_it_t map_it_t;
    typedef std::map<std::pair<mpo_index_t, mpo_index_t>, block_matrix<Matrix, SymmGroup>, compressor_detail::pair_cmp_inv> map_t;
    typedef std::vector<std::pair<mpo_index_t, typename Matrix::value_type> > ret_t;

public:
    column_iterator(MPOTensor<Matrix, SymmGroup> const & mpos, charge lc,
                    Index<SymmGroup> const & ii, bond_charge_t const & bd, Index<SymmGroup> const & phys, charge_overlap<Matrix, SymmGroup> co)
        : cur_col(-1), c(-1), coldim(0), mpo_source(mpos), left_charge(lc), in_index(ii), out_bond_charge(bd), phys_i(phys), op_offset(co)
    {
        // variable init
        //empty_count = 0; // debug
        block_pos = 0; block_size = 0;

        // create set with physical charges
        for (int k = 0; k < phys_i.size(); ++k)
            phys_set.insert(phys_i[k].first);

        // column init
        for (mpo_index_t c = 0; c < mpo_source.col_dim(); ++c) {  
            coldim += op_offset(left_charge, out_bond_charge[c]);
        }

        // row init
        rowdim = in_index.size_of_block(left_charge);
        in_index[in_index.position(SymmGroup::IdentityCharge)].second -= 2;
        rowstart = 2;
        for (typename Index<SymmGroup>::iterator it = in_index.begin(); it != in_index.end(); ++it)
        {
            if (it->first > left_charge)
                rowstart += it->second;
            else {
                rowend = rowstart + it->second;
                break;
            }
        }

        // resorted (local) MPOTensor for column-major-like lookup
        if (left_charge == SymmGroup::IdentityCharge) {

            for (map_it_t it = mpo_source.data_.lower_bound(std::make_pair(0, 0));
                it != mpo_source.data_.lower_bound(std::make_pair(2,0)); ++it)
                    col_order[it->first] = it->second;

            for (map_it_t it = mpo_source.data_.lower_bound(std::make_pair(rowstart, 0));
                it != mpo_source.data_.lower_bound(std::make_pair(rowend+1,0)); ++it) {
                    typename map_t::key_type key_ = std::make_pair(it->first.first-rowstart+2, it->first.second);
                    col_order[key_] = it->second;
            }
        }
        else {

            for (map_it_t it = mpo_source.data_.lower_bound(std::make_pair(rowstart, 0));
                it != mpo_source.data_.lower_bound(std::make_pair(rowend+1,0)); ++it) {
                    typename map_t::key_type key_ = std::make_pair(it->first.first-rowstart, it->first.second);
                    col_order[key_] = it->second;
            }
        }

        ++(*this); // advance to first valid position
    }

    column_iterator & operator++()
    {
        // operator block
        block_pos++;
        if (block_pos >= block_size) {
            do {
                c++;
                block_size = op_offset(left_charge, out_bond_charge[c]);
            } while(block_size == 0);
            block_pos = 0;
            op_it_low = col_order.lower_bound(std::make_pair(0,c));
            op_it_high = col_order.lower_bound(std::make_pair(0,c+1));
        }
        op_it_loc = op_it_low;

        // intra operator block position (0,1,2,3)
        charge delta = SymmGroup::fuse(left_charge, -out_bond_charge[c]); 
        local_ls = 0;
        for (mpo_index_t ls = 0; ls < phys_i.size(); ++ls) {
            op_out = SymmGroup::fuse(phys_i[ls].first, delta);
            if (phys_set.count(op_out) == 0)
                continue;

            if (local_ls == block_pos) {
                local_ls = ls; 
                break;
            }
            local_ls++;
        }

        cur_col++;
        return *this;
    }

    ret_t operator*() const 
    {
        ret_t ret;

        while(op_it_loc != op_it_high)
        {
            mpo_index_t rloc = op_it_loc->first.first;
            block_matrix<Matrix, SymmGroup> const & op = op_it_loc->second;

            if ( op.has_block(phys_i[local_ls].first, op_out) ) {
                mpo_index_t cs = op.left_basis().position(phys_i[local_ls].first);
                ret.push_back(std::make_pair(rloc, op[cs](0,0)));
            }

            //mpo_index_t r = ( (left_charge == SymmGroup::IdentityCharge) && (rloc < 2) ? rloc : rloc+rowstart);
            //maquis::cout << left_charge << "(" << rloc << "," << cur_col << "): "
            //             << ret[rloc] << " <- (" << r << "," << c << ")" << std::endl;
            ++op_it_loc;
        }

        // Debug
        //if (!has_op) {
            //maquis::cout << "has no op! c=" << c << " block_size: "
            //             << op_offset(left_charge, out_bond_charge[c]) << std::endl;

            //for (typename MPOTensor<Matrix, SymmGroup>::data_t::const_iterator it = mpo_source.data_.begin();
            //              it != mpo_source.data_.end(); ++it)
            //    if( it->first.second == c)
            //        maquis::cout << "Note: rrange: " << rowstart << "-" << rowend <<
            //        ", has (" << it->first.first << "," << c << "): " << it->second << std::endl;
        //    empty_count++;
        //}

        return ret;
    }

    mpo_index_t key() const {
        return (mpo_index_t(c) << 4) + (local_ls << 2) + mpo_index_t(phys_i.position(op_out));
    }

    inline bool end() const { return cur_col == coldim; }
    inline mpo_index_t col() const { return cur_col; }
    inline mpo_index_t col_dim() const { return coldim; }
    inline mpo_index_t row_dim() const { return rowdim; }
    //inline mpo_index_t empty_cols() const { return empty_count; }

private:
    // Iteration state
    offset_t cur_col, c;
    mpo_index_t block_pos, block_size, local_ls;
    charge op_out;
    typename map_t::const_iterator op_it_low, op_it_high;
    mutable typename map_t::const_iterator op_it_loc;

    // Constant state data
    mpo_index_t rowstart, rowend, rowdim, coldim;
    map_t col_order;
    std::set<charge> phys_set;

    // Input data
    charge left_charge;
    MPOTensor<Matrix, SymmGroup> const & mpo_source;
    Index<SymmGroup> phys_i;
    Index<SymmGroup> in_index;
    bond_charge_t const & out_bond_charge;
    mutable charge_overlap<Matrix, SymmGroup> op_offset;

    // Debug
    //mutable mpo_index_t empty_count;
};

///////////////////////////////////////////////////////////////////////////////

template<class Matrix, class SymmGroup>
class compressor
{
    typedef typename SymmGroup::charge charge;
    typedef typename alps::numeric::associated_dense_matrix<Matrix>::type dense_matrix;
    typedef typename MPO<Matrix, SymmGroup>::elem_type elem_type;
    typedef typename MPOIndexer<Matrix, SymmGroup>::map_it_t map_it_t;
    typedef typename MPOIndexer<Matrix, SymmGroup>::mpo_index_t mpo_index_t;

public:
    compressor() {}
    compressor(Index<SymmGroup> const & phys) : phys_i(phys), op_offset(phys_i)
    {
        charge_deltas = op_offset.get_charge_deltas();
        for (int ii=0; ii<phys_i.size(); ++ii)
            phys_set.insert(phys_i[ii].first);
    }

    void compress(MPO<Matrix, SymmGroup> const & mpo_in,
                  MPO<Matrix, SymmGroup> & mpo_out, double cutoff);

    template <class Matrix2>
    static void convert_to_dense_matrix(MPO<Matrix, SymmGroup> const & mpo_in, MPO<Matrix2, SymmGroup> & mpo_out);

private:
    block_matrix<dense_matrix, SymmGroup> make_left_matrix(MPO<Matrix, SymmGroup> const & mpo_in, std::size_t p);
    block_matrix<dense_matrix, SymmGroup> make_right_matrix(MPO<Matrix, SymmGroup> const & mpo_in, std::size_t p);
    block_matrix<dense_matrix, SymmGroup> make_M_matrix(MPO<Matrix, SymmGroup> const & mpo_in,
                                                        MPO<Matrix, SymmGroup> & mpo_out, std::size_t p);

    void replace_pair(block_matrix<dense_matrix, SymmGroup> const & left,
                      block_matrix<dense_matrix, SymmGroup> const & right,
                      MPO<Matrix, SymmGroup> const & mpo_in, MPO<Matrix, SymmGroup> & mpo_out, std::size_t p);

    Index<SymmGroup> phys_i;
    std::set<charge> phys_set;

    std::vector<charge> charge_deltas;
    std::map<charge, std::vector<mpo_index_t> > M_keys;
    
    charge_overlap<Matrix, SymmGroup> op_offset;
    MPOIndexer<Matrix, SymmGroup> bond_descriptor;
};

template <class Matrix, class SymmGroup>
void compressor<Matrix, SymmGroup>::compress(MPO<Matrix, SymmGroup> const & mpo_in, MPO<Matrix, SymmGroup> & mpo_out, double cutoff)
{
    bond_descriptor = MPOIndexer<Matrix, SymmGroup>(mpo_in);
    mpo_out.resize(mpo_in.size());

    ////////////////////////////////////////////////////////////////
    // Playground
    /*
    unsigned probe = 3;
    column_iterator<Matrix, SymmGroup> p1adaptor(mpo_in[probe], bond_descriptor.index(probe).charges()[6], 
                                                bond_descriptor.index(probe), bond_descriptor[probe+1], phys_i, op_offset);

    //maquis::cout << bond_descriptor.index(probe+1);
    
    for (int cc = 0; cc < p1adaptor.col_dim(); ++cc)
    {
        std::vector<typename Matrix::value_type> bla = *p1adaptor;
        ++p1adaptor;
    }
    maquis::cout << "Matrix block " << p1adaptor.row_dim() << "x" << p1adaptor.col_dim() << ", empty cols: " << p1adaptor.empty_cols() << std::endl;
    */
    ////////////////////////////////////////////////////////////////

    std::cout << "starting MPO SVD compression\n";
    
    Timer cot("compression"), lm("left_matrix"), rm("right_matrix"), mm("m_matrix"), tsvd("SVD"), rep("replace");
    cot.begin();

    for (int p = 0; p < mpo_in.size()-1; ++p) {

        //maquis::cout << "Stalling after replace @ p=" << p << std::endl;
        //usleep(7*1000000);
        
        mm.begin();
        block_matrix<dense_matrix, SymmGroup> MM = make_M_matrix(mpo_in, mpo_out, p);
        mm.end();

        //maquis::cout << "Stalling after M @ p=" << p << std::endl;
        //usleep(7*1000000);

        block_matrix<dense_matrix, SymmGroup> U, V;
        block_matrix<typename alps::numeric::associated_real_diagonal_matrix<dense_matrix>::type, SymmGroup> S, Sqrt;

        //maquis::cout << "MM shape @ " << p << std::endl; mem_footprint(MM); maquis::cout << std::endl;
        tsvd.begin();
        maquis::cout << "starting svd_truncate @ " << p << std::endl;
        svd_truncate(MM, U, V, S, cutoff, 100000, false);
        tsvd.end();

        Sqrt = sqrt(S);
        for (mpo_index_t k = 0; k < Sqrt.n_blocks(); ++k) {
            for (mpo_index_t s = 0; s < Sqrt[k].num_cols(); ++s)
                for( mpo_index_t j = 0; j < U[k].num_rows(); ++j)
                    U[k](j,s) *= Sqrt[k](s,s);

            for( mpo_index_t j = 0; j < V[k].num_cols(); ++j)
                for (mpo_index_t s = 0; s < Sqrt[k].num_cols(); ++s)
                    V[k](s,j) *= Sqrt[k](s,s);
        }
        //gemm(U, Sqrt, left);
        //gemm(Sqrt, V, right);

        //maquis::cout << "Stalling after SVD t @ p=" << p << std::endl;
        //usleep(7*1000000);
        
        maquis::cout << "MPO bond truncation: " << bond_descriptor.index(p+1).sum_of_sizes() << " -> ";
        rep.begin();
        replace_pair(U, V, mpo_in, mpo_out, p);
        rep.end();
        maquis::cout << bond_descriptor.index(p+1).sum_of_sizes() << std::endl;
        maquis::cout << std::endl;
    }
    cot.end();
    maquis::cout << "Rest: " << cot.get_time() - tsvd.get_time() - rep.get_time() - mm.get_time() << std::endl;
}

template <class Matrix, class SymmGroup>
block_matrix<typename compressor<Matrix, SymmGroup>::dense_matrix, SymmGroup> compressor<Matrix, SymmGroup>::make_left_matrix(MPO<Matrix, SymmGroup> const & mpo_in, std::size_t p)
{
    Index<SymmGroup> left_i = phys_i * adjoin(phys_i) * bond_descriptor.index(p);
    Index<SymmGroup> right_i = bond_descriptor.index(p+1);
    left_i = common_subset(left_i, right_i);
    block_matrix<dense_matrix, SymmGroup> ret(left_i, right_i);

    std::map<charge, std::vector<int> > outr_lookup;
    for (int k = 0; k < charge_deltas.size(); ++k) {
        int outr = -1;
        charge rc = charge_deltas[k];
        std::vector<int> & vec_ref = outr_lookup[rc];

        for (int r = 0; r < mpo_in[p].row_dim(); ++r)
        {  
            vec_ref.push_back(outr);
            outr += op_offset(bond_descriptor[p][r], rc);
        }
    }

    typedef std::set<std::pair<mpo_index_t, mpo_index_t>, compressor_detail::pair_cmp_inv> rset_t;
    rset_t col_order;
    for (map_it_t it = mpo_in[p].data_.begin(); it != mpo_in[p].data_.end(); ++it)
        col_order.insert(it->first);

    std::map<charge, mpo_index_t> visited_c_basis;
    charge rc = bond_descriptor[p+1][0], rcprev = rc;
    mpo_index_t cprev = 0;
    for(typename rset_t::iterator it=col_order.begin(); it != col_order.end(); ++it) {
        mpo_index_t r = it->first, c = it->second;

        if (c != cprev) {
            rc = bond_descriptor[p+1][c];
            visited_c_basis[rcprev]++;
            rcprev = rc;
        }
        cprev = c; 
        int outr = outr_lookup[rc][r];
        block_matrix<Matrix, SymmGroup> const & op = mpo_in[p](r,c);

        charge delta = SymmGroup::fuse(bond_descriptor[p][r], -rc);
        for (mpo_index_t ls = 0; ls < phys_i.size(); ++ls) {

                charge op_out = SymmGroup::fuse(phys_i[ls].first, delta);
                if (phys_set.count(op_out) == 0)
                    continue;
                
                outr++;
                
                if (! op.has_block(phys_i[ls].first, op_out) )
                    continue;                       
                
                mpo_index_t cs = op.left_basis().position(phys_i[ls].first);
                ret(std::make_pair(rc, outr),
                    std::make_pair(rc, visited_c_basis[rc])) = op[cs](0,0);
                
        }
    }
    
    return ret;
}

/*  Superseded by column_iterator (direct generation in make_M_matrix)
template <class Matrix, class SymmGroup>
block_matrix<typename compressor<Matrix, SymmGroup>::dense_matrix, SymmGroup> compressor<Matrix, SymmGroup>::make_right_matrix(MPO<Matrix, SymmGroup> const & mpo_in, std::size_t p)
{
    Index<SymmGroup> left_i = bond_descriptor.index(p);
    Index<SymmGroup> right_i = adjoin(phys_i) * phys_i * bond_descriptor.index(p+1);
    right_i = common_subset(right_i, left_i);
    block_matrix<dense_matrix, SymmGroup> ret(left_i, right_i);

    std::map<charge, std::vector<int> > outc_lookup;
    for (int k = 0; k < charge_deltas.size(); ++k) {
        int outc = -1;
        charge lc = charge_deltas[k];
        std::vector<int> & vec_ref = outc_lookup[lc];

        for (int c = 0; c < mpo_in[p].col_dim(); ++c)
        {  
            vec_ref.push_back(outc);
            outc += op_offset(bond_descriptor[p+1][c], lc);
        }
    }


    std::map<charge, size_t> visited_r_basis;
    charge lc = bond_descriptor[p][0], lcprev = lc;
    std::size_t rprev = 0;
    for(map_it_t it=mpo_in[p].data_.begin(); it != mpo_in[p].data_.end(); ++it) {
        std::size_t r = it->first.first, c = it->first.second;
        
        if (r != rprev) {
            lc = bond_descriptor[p][r];
            visited_r_basis[lcprev]++;
            lcprev = lc;
        }
        rprev = r;
        int outc = outc_lookup[lc][c];
        block_matrix<Matrix, SymmGroup> const & op = it->second;
        
        charge delta = SymmGroup::fuse(lc, -bond_descriptor[p+1][c]);
        for (size_t ls = 0; ls < phys_i.size(); ++ls) {

                charge op_out = SymmGroup::fuse(phys_i[ls].first, delta);
                if (phys_set.count(op_out) == 0)
                    continue;
                
                outc++;
                
                if (! op.has_block(phys_i[ls].first, op_out) )
                    continue;
                
                std::size_t cs = op.left_basis().position(phys_i[ls].first);
                ret(std::make_pair(lc, visited_r_basis[lc]),
                    std::make_pair(lc, outc)) = op[cs](0,0);
                
        }
    }
    
    return ret;
}
*/

template <class Value_Type>
inline void axpy(const int*, const Value_Type*, const Value_Type*, const int*, Value_Type*, const int*)
{ throw std::runtime_error("Blas axpy call not implemented for Matrix::value_type used\n"); }

template < >
inline void axpy(const int* nr, const double* val, const double* x, const int* ix, double* y, const int* iy)
{ daxpy_(nr, val, x, ix, y, iy); }

template < >
inline void axpy(const int* nr, const std::complex<double>* val, const std::complex<double>* x, const int* ix, std::complex<double>* y, const int* iy)
{ zaxpy_(nr, val, x, ix, y, iy); }

template <class Matrix, class SymmGroup> block_matrix<typename compressor<Matrix, SymmGroup>::dense_matrix, SymmGroup>
compressor<Matrix, SymmGroup>::make_M_matrix(MPO<Matrix, SymmGroup> const & mpo_in, MPO<Matrix, SymmGroup> & mpo_out, std::size_t p)
{
    Index<SymmGroup> left_i = phys_i * adjoin(phys_i) * bond_descriptor.index(p);
    Index<SymmGroup> mid_i1 = bond_descriptor.index(p+1);
    Index<SymmGroup> mid_i2 = bond_descriptor.index(p+1);
    Index<SymmGroup> right_i = adjoin(phys_i) * phys_i * bond_descriptor.index(p+2);

    left_i = common_subset(left_i, mid_i1);
    right_i = common_subset(right_i, mid_i2);

    if (left_i.size() != right_i.size()) { maquis::cout << "M index left != right in size\n"; }
    block_matrix<dense_matrix, SymmGroup> left = make_left_matrix( (p>0) ? mpo_out : mpo_in, p);
    //maquis::cout << "L shape @ " << p << std::endl; mem_footprint(left); maquis::cout << std::endl;

    block_matrix<dense_matrix, SymmGroup> ret;

    // charge loop, 13 main blocks in chem
    for (mpo_index_t csl = 0; csl < left_i.size(); ++csl)
    {
        std::vector<mpo_index_t> & locKey= M_keys[left_i[csl].first];
        locKey.clear();

        column_iterator<Matrix, SymmGroup> rAdaptor(mpo_in[p+1], left_i[csl].first,
                                       bond_descriptor.index(p+1), bond_descriptor[p+2], phys_i, op_offset);

        //dense_matrix block(left_i[csl].second, right_i[csl].second);
        mpo_index_t newpos = ret.insert_block(dense_matrix(left_i[csl].second, right_i[csl].second),
                                              left_i[csl].first, left_i[csl].first);
        dense_matrix & block = ret[newpos];

        // right columns
        mpo_index_t colcnt = 0;
        for( ; !rAdaptor.end(); ++rAdaptor)
        {
            // left columns / right row
            //int nr = left[csl].num_rows(), one=1;
            //mpo_index_t k = rAdaptor.col();
            std::vector<std::pair<mpo_index_t, typename Matrix::value_type> > mult = (*rAdaptor);
            for (typename std::vector<std::pair<mpo_index_t, typename Matrix::value_type> >::const_iterator col_it
                    = mult.begin(); col_it != mult.end(); ++col_it)
            {
                mpo_index_t s = col_it->first;
                typename Matrix::value_type val = col_it->second;

                // left row
                for (mpo_index_t l=0; l < left[csl].num_rows(); ++l)
                    block(l, colcnt) += left[csl](l,s) * val;

                //axpy(&nr, &val, left[csl].col(s).first, &one, ret[csl].col(colcnt).first, &one);
            }

            if (mult.size() > 0) {
                colcnt++;
                locKey.push_back(rAdaptor.key());
            }

            //mpo_index_t key = rAdaptor.key();
            //maquis::cout << "site " << p << "k" << colcnt << "/" << k <<
            //    " key: " << (key >> 4) << "," << ((key&12u) >> 2) << "," << (key&3u) << std::endl;
        }
        //maquis::cout << "resized " << p << "@" << left_i[csl].first << " block " << right_i[csl].second << " -> " << colcnt << std::endl;
        ret.resize_block(left_i[csl].first, left_i[csl].first, left_i[csl].second, colcnt, false);

        /*
        maquis::cout << "row profiles: ";
        for (int r=0; r < block.num_rows(); ++r) {
            typedef typename dense_matrix::const_row_element_iterator cri;
            std::pair<cri, cri> pr = block.row(r);
            mpo_index_t non_zero_cnt=0;
            for(;pr.first != pr.second; ++pr.first)
                if ( std::abs(*pr.first) > 1e-40 ) non_zero_cnt++;

            maquis::cout << non_zero_cnt << " ";
        }
        maquis::cout << std::endl;
        */

        //maquis::cout << "block num_cols " << block.num_cols() << " is shrinkable " << block.is_shrinkable() << std::endl;

        // M-block built

        //dense_matrix U, V, sU, sV;
        //typename alps::numeric::associated_real_diagonal_matrix<dense_matrix>::type S;

        //svd(block, U, V, S, 1e-12, 100000, false);
        //sqrt(S);

        //sU.resize(U.num_rows(), U.num_cols());
        //sV.resize(V.num_rows(), V.num_cols());
        //gemm(U, S, sU);
        //gemm(S, V, sV);
        //
        //replace_pair(sU, sV, mpo_in, mpo_out, p);
    } 
    
    return ret;
}

template <class Matrix, class SymmGroup>
void compressor<Matrix, SymmGroup>::replace_pair(block_matrix<typename compressor<Matrix, SymmGroup>::dense_matrix, SymmGroup> const & left,
                                             block_matrix<typename compressor<Matrix, SymmGroup>::dense_matrix, SymmGroup> const & right,
                                             MPO<Matrix, SymmGroup> const & mpo_in,
                                             MPO<Matrix, SymmGroup> & mpo_out,
                                             std::size_t p)
{
    assert( left.right_basis() == right.left_basis() );
    bond_descriptor.update_bond(left.right_basis(), p+1);
    
    mpo_out[p] = MPOTensor<Matrix, SymmGroup>( (p>0) ? mpo_out[p].row_dim() : mpo_in[p].row_dim(),
                                              bond_descriptor.index(p+1).sum_of_sizes());
    
    mpo_index_t term_cnt=0;
    std::size_t bstart = 0, bend = 0;
    for ( typename Index<SymmGroup>::const_iterator it = left.right_basis().begin();
        it != left.right_basis().end(); ++it) {
        bend += it->second; 
        charge rc = bond_descriptor[p+1][bstart];
        dense_matrix const & block = left(rc, rc);

        std::size_t cloc = 0;
        for (size_t c = bstart; c < bend; ++c) {
            int outr = -1;
            for (size_t r = 0; r < mpo_out[p].row_dim(); ++r) {

                charge delta = SymmGroup::fuse(bond_descriptor[p][r], -rc);
                for (size_t ls = 0; ls < phys_i.size(); ++ls)
                {
                    charge op_out = SymmGroup::fuse(phys_i[ls].first, delta);
                    if (phys_set.count(op_out) == 0)
                        continue;
                    
                    outr++;
                    
                    typename Matrix::value_type val = block(outr, cloc);
                    if (std::abs(val) > 1e-40) {
                        term_cnt++;
                        mpo_out[p](r,c).insert_block(Matrix(1,1,val), phys_i[ls].first, op_out);
                    }
                }
            }
            cloc++;
        }
        bstart = bend; 
    }
    //maquis::cout << "left replace @" << p << ": " << term_cnt << " terms in " << mpo_out[p].data_.size() << " map blocks\n";
    //maquis::cout << "original size " << mpo_in[p].data_.size() << " map blocks\n";
    
    mpo_out[p+1] = MPOTensor<Matrix, SymmGroup>(bond_descriptor.index(p+1).sum_of_sizes(),
                                                mpo_in[p+1].col_dim());
    
    term_cnt=0;
    bstart = 0, bend = 0;
    for ( typename Index<SymmGroup>::const_iterator it = left.right_basis().begin();
        it != left.right_basis().end(); ++it) {
        bend += it->second; 
        charge lc = bond_descriptor[p+1][bstart];
        dense_matrix const & block = right(lc, lc);
        std::vector<mpo_index_t> const & locKey = M_keys[lc];

        std::size_t rloc = 0;
        for (size_t r = bstart; r < bend; ++r) {
            for (size_t k = 0; k < locKey.size(); ++k) {

                typename Matrix::value_type val = block(rloc, k);
                if (std::abs(val) > 1e-40) {
                    term_cnt++;
                    mpo_index_t key = locKey[k];
                    mpo_index_t c = (key >> 4);
                    mpo_out[p+1](r,c).insert_block(Matrix(1,1,val), phys_i[((key&12u)>>2)].first, phys_i[key&3u].first);
                    // ***** Debug *********
                    //maquis::cout << "right checkback " << p+1 << "(" << r << "," << c << ")[" << ls << "," << phys_i.position(op_out) << "]: " << val
                    //             << " | key " << (key >> 4) << "," << ((key&12u)>>2) << "," << (key&3u) << std::endl;
                }
            }
            rloc++;
        }
        bstart = bend; 
    }
    //maquis::cout << "right replace @" << (p+1) << ": " << term_cnt << " terms in " << mpo_out[p+1].data_.size() << " map blocks\n";
    //maquis::cout << "original size " << mpo_in[p+1].data_.size() << " map blocks\n";

}

template <class Matrix, class SymmGroup>
template <class Matrix2>
void compressor<Matrix, SymmGroup>::convert_to_dense_matrix(MPO<Matrix, SymmGroup> const & mpo_in, MPO<Matrix2, SymmGroup> & mpo_out)
{
    mpo_out.resize(mpo_in.size());
    for (int p = 0; p < mpo_in.size(); ++p)
    {
        mpo_out[p] = MPOTensor<Matrix2, SymmGroup>(mpo_in[p].row_dim(), mpo_in[p].col_dim());
        for (map_it_t it = mpo_in[p].data_.begin(); it != mpo_in[p].data_.end(); ++it)
        {
            block_matrix<Matrix2, SymmGroup> new_op;
            block_matrix<Matrix, SymmGroup> const & old_op = it->second;
            for(int b=0; b < old_op.n_blocks(); ++b)
                new_op.insert_block( Matrix2(1,1, old_op[b](0,0)),
                                     old_op.left_basis()[b].first, old_op.right_basis()[b].first);

            mpo_out[p](it->first.first, it->first.second) = new_op;
        } 
    }
}

// Relabel the MPO-indices such that the 13 possible bond charges are grouped together
template <class Tuple, class Matrix, class SymmGroup> void
charge_sort(std::vector< std::vector<
                //boost::tuple<std::size_t, std::size_t, block_matrix<Matrix, SymmGroup> > > > & prempo,
                Tuple > > & prempo,
            MPO<Matrix, SymmGroup> const & mpo_in)
{
    //typedef boost::tuple<std::size_t, std::size_t, block_matrix<Matrix, SymmGroup> > block;
    typedef typename SymmGroup::charge charge;

    Timer csort("MPO_charge_sort"); csort.begin();
    MPOIndexer<Matrix, SymmGroup> mpo_in_index(mpo_in);
    for (std::size_t p = 0; p < prempo.size()-1; ++p)
    {
        // TODO: Check types
        std::map<std::size_t, std::size_t> reorder_map;
        std::map<charge, std::size_t> visited_delta_basis;
        std::map<charge, std::size_t> bond_index_offsets;

        // determine reordering
        for (typename std::vector<Tuple>::iterator it = prempo[p].begin();
                it != prempo[p].end(); ++it)
        {
            std::size_t c = boost::tuples::get<1>(*it);
            if (c < 2) continue;
            charge out_charge = mpo_in_index[p+1][c];
            reorder_map[c] = visited_delta_basis[out_charge]++;
        }
        
        // calculate accumulated charge sector sizes
        std::map<charge, std::size_t> offsets;
        std::size_t offset = 2;
        for (typename std::map<charge, std::size_t>::reverse_iterator it = visited_delta_basis.rbegin();
                it != visited_delta_basis.rend(); ++it)
        {
            offsets[it->first] = offset;
            offset += it->second; 
        }

        // add accumulated offsets to the reordering
        for (typename std::vector<Tuple>::iterator it = prempo[p].begin();
                it != prempo[p].end(); ++it)
        {
            std::size_t c = boost::tuples::get<1>(*it);
            if (c < 2) continue;
            charge out_charge = mpo_in_index[p+1][c];
            reorder_map[c] += offsets[out_charge];
        }

        // carry out reordering
        for (typename std::vector<Tuple>::iterator it = prempo[p].begin(); it != prempo[p].end(); ++it)
            if (reorder_map.count(boost::tuples::get<1>(*it)) > 0)
                boost::tuples::get<1>(*it) = reorder_map[boost::tuples::get<1>(*it)];
        for (typename std::vector<Tuple>::iterator it = prempo[p+1].begin(); it != prempo[p+1].end(); ++it)
            if (reorder_map.count(boost::tuples::get<0>(*it)) > 0)
                boost::tuples::get<0>(*it) = reorder_map[boost::tuples::get<0>(*it)];

    }
    csort.end();
}


#endif
