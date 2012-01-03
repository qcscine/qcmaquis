
#ifndef MPO_MANIP_H
#define MPO_MANIP_H

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/multi_index.h"

#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/mp_tensors/mpo.h"


template <class Matrix, class SymmGroup>
MPO<Matrix, SymmGroup> block_to_mpo(Index<SymmGroup> const & phys_i,
                                    block_matrix<Matrix, SymmGroup> block,
                                    std::size_t length)
{
    MPO<Matrix, SymmGroup> mpo(length);
    
    typedef typename MultiIndex<SymmGroup>::index_id index_id;
    typedef typename MultiIndex<SymmGroup>::set_id set_id;
    
    
    Index<SymmGroup> alpha_i, beta_i;
    alpha_i.insert( std::make_pair(SymmGroup::IdentityCharge,1) );
    beta_i.insert( std::make_pair(SymmGroup::IdentityCharge,1) );
    
    
    MultiIndex<SymmGroup> midx;
    index_id alpha, beta;
    std::vector<std::pair<index_id, bool> > in_left, in_right;
    std::vector<std::pair<index_id, index_id> > phys_ids;
    
    alpha = midx.insert_index(alpha_i);
    for (size_t i=0; i<length; ++i) {
        index_id id1 = midx.insert_index(phys_i);
        index_id id2 = midx.insert_index(phys_i);
        phys_ids.push_back( std::make_pair(id1, id2) );
        
        in_left.push_back( std::make_pair(id1, true) );
        in_right.push_back( std::make_pair(id2, true) );
    }
    beta = midx.insert_index(beta_i);
    in_left.push_back( std::make_pair(beta, true) );
    in_right.push_back( std::make_pair(alpha, true) );
    
    set_id curr_s = midx.create_set(in_left, in_right);

    
    for (size_t p=0; p<length; ++p) {
        assert( phys_ids.size() == length-p );
        

        block_matrix<Matrix, SymmGroup> btmp;
        
        // Reshape & SVD
        if (p < length-1) {
            std::vector<std::pair<index_id, bool> > svd_left, svd_right;
            
            for (size_t i=0; i<length-p; ++i) {
                if (i == 0) {
                    svd_left.push_back( std::make_pair(phys_ids[i].first, true) );
                    svd_left.push_back( std::make_pair(phys_ids[i].second, false) );
                } else {
                    svd_right.push_back( std::make_pair(phys_ids[i].first, false) );
                    svd_right.push_back( std::make_pair(phys_ids[i].second, true) );
                }
            }
            svd_left.push_back( std::make_pair(alpha, false) );
            svd_right.push_back( std::make_pair(beta, false) );
            
            set_id svd_s = midx.create_set(svd_left, svd_right);
            
            reshape(midx, curr_s, svd_s, block, btmp);
            
            block_matrix<Matrix, SymmGroup> U, V;
            block_matrix<typename maquis::types::associated_diagonal_matrix<Matrix>::type, SymmGroup> S, Ssqrt;
            svd(btmp, U, V, S);
            Ssqrt = sqrt(S);
            gemm(U, Ssqrt, btmp);
            gemm(Ssqrt, V, block);
        } else {
            assert( phys_ids.size() == 1 );
            std::vector<std::pair<index_id, bool> > vec_left, vec_right;
            index_id id1 = phys_ids[0].first;
            index_id id2 = phys_ids[0].second;
            vec_left.push_back( std::make_pair(id1, true) );
            vec_left.push_back( std::make_pair(id2, false) );
            vec_left.push_back( std::make_pair(alpha, false) );
            vec_right.push_back( std::make_pair(beta, false) );
            
            set_id new_s = midx.create_set(vec_left, vec_right);

            reshape(midx, curr_s, new_s, block, btmp);
        }
        
        // Reshaping btmp in MPOTensor
        Index<SymmGroup> const & aux_left_i = alpha_i;
        Index<SymmGroup> const & aux_right_i = adjoin(btmp.right_basis());
        MultiIndex<SymmGroup> midx_mpo;
        index_id aux_left = midx_mpo.insert_index(aux_left_i);
        index_id phys_mpo1 = midx_mpo.insert_index(phys_i);
        index_id phys_mpo2 = midx_mpo.insert_index(phys_i);
        index_id aux_right = midx_mpo.insert_index(aux_right_i);
        
        std::vector<std::pair<index_id, bool> > vec_curr_mpo_l, vec_curr_mpo_r;
        vec_curr_mpo_l.push_back( std::make_pair(phys_mpo1, true) );
        vec_curr_mpo_l.push_back( std::make_pair(phys_mpo2, false) );
        vec_curr_mpo_l.push_back( std::make_pair(aux_left, false) );
        vec_curr_mpo_r.push_back( std::make_pair(aux_right, false) );
        
        set_id curr_mpo_s = midx_mpo.create_set(vec_curr_mpo_l, vec_curr_mpo_r);
        mpo[p] = MPOTensor<Matrix, SymmGroup>(aux_left_i.sum_of_sizes(), aux_right_i.sum_of_sizes());
        for (short run=0; run<2; ++run) {
            
            if (run == 1)
                for (size_t r=0; r<mpo[p].row_dim(); ++r)
                    for (size_t c=0; c<mpo[p].col_dim(); ++c)
                        if (mpo[p].has(r, c))
                            mpo[p](r, c).allocate_blocks();
            
            for (index_product_iterator<SymmGroup> it = midx_mpo.begin();
                 it != midx_mpo.end();
                 it++)
            {
                typename MultiIndex<SymmGroup>::coord_t lc, rc;
                boost::tie(lc, rc) = midx_mpo.get_coords(curr_mpo_s, *it);
                if (!btmp.has_block(lc.first, rc.first))
                    continue;

                size_t li = aux_left_i.position( (*it)[aux_left] );
                size_t ri = aux_right_i.position( (*it)[aux_right] );
                
                if (btmp(lc, rc) != 0.) {
                    if (run == 0)
                        mpo[p](li, ri).reserve((*it)[phys_mpo1].first, (*it)[phys_mpo2].first,
                                               phys_i.size_of_block((*it)[phys_mpo1].first), phys_i.size_of_block((*it)[phys_mpo2].first));
                    else
                        mpo[p](li, ri)((*it)[phys_mpo1], (*it)[phys_mpo2]) = btmp(lc, rc);
                }
            }
        }
        
        
        // Preparing new loop
        if (p < length-1) {
            midx.clear();
            phys_ids.clear();
            std::vector<std::pair<index_id, bool> > out_left, out_right;
            
            alpha_i = adjoin(block.left_basis());
            alpha = midx.insert_index(alpha_i);
            for (size_t i=p+1; i<length; ++i) {
                index_id id1 = midx.insert_index(phys_i);
                index_id id2 = midx.insert_index(phys_i);
                phys_ids.push_back( std::make_pair(id1, id2) );
                
                out_right.push_back( std::make_pair(id1, false) );
                out_right.push_back( std::make_pair(id2, true) );
            }
            beta = midx.insert_index(beta_i);
            out_left.push_back( std::make_pair(alpha, false) );
            out_right.push_back( std::make_pair(beta, false) );
            
            curr_s = midx.create_set(out_left, out_right);
        }
    }
    
    return mpo;
}


#endif
