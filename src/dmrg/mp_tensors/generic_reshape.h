
#ifndef GENERIC_RESHAPE_H
#define GENERIC_RESHAPE_H

#include "block_matrix/block_matrix.h"
#include "block_matrix/multi_index.h"


template<class Matrix, class SymmGroup>
void reshape(MultiIndex<SymmGroup> const & midx,
             typename MultiIndex<SymmGroup>::set_id in_set,
             typename MultiIndex<SymmGroup>::set_id out_set,
             block_matrix<Matrix, SymmGroup> const & m1,
             block_matrix<Matrix, SymmGroup> & m2)
{   
    static Timer timer("generic reshape");
    timer.begin();
    
    
    m2 = block_matrix<Matrix, SymmGroup>();
    
    typedef std::size_t size_t;
    typedef typename SymmGroup::charge charge;
    typedef typename MultiIndex<SymmGroup>::coord_t coord_t;
    
    for (int run = 0; run < 2; ++run) {
        bool pretend = (run == 0);

//        std::cout << "pretend == " << pretend << std::endl;
        
        if (!pretend)
            m2.allocate_blocks();

//        if (!pretend)
//            std::cout << "m2 after allocating: " << m2.description() << std::endl;

        
        for (size_t block = 0; block < m1.n_blocks(); ++block)
        {
//            std::cout << "block " << block << " -- start" <<std::endl;
            
            charge in_l_charge = m1.left_basis()[block].first;
            charge in_r_charge = m1.right_basis()[block].first;
            
//            std::cout << "block " << block << " -- sector(" << in_l_charge << "," << in_r_charge << ")" <<std::endl;
            
            Matrix const & in_block = m1[block];
//            Matrix const & in_block = m1(in_l_charge, in_r_charge);
            
            for (size_t i=0; i<num_rows(in_block); ++i)
                for (size_t j=0; j<num_cols(in_block); ++j) {
                    coord_t in_left = std::make_pair(in_l_charge, i);
                    coord_t in_right = std::make_pair(in_r_charge, j);
                    coord_t out_left, out_right;
//                    std::cout << "block[" << block << "](" << i << "," << j << ") -- do convert" <<std::endl;
                    boost::tie(out_left, out_right) = midx.convert_coords(in_set, in_left, in_right, out_set);
                    
//                    std::cout << "block[" << block << "](" << i << "," << j << ") -- out_left " << out_left << std::endl;
//                    std::cout << "block[" << block << "](" << i << "," << j << ") -- out_right " << out_right << std::endl;
//                    std::cout << "block[" << block << "](" << i << "," << j << ") -- do reserve/assign" <<std::endl;
                    if (pretend)
                        m2.reserve_pos(out_left.first, out_right.first,
                                       out_left.second, out_right.second);
                    else
                        m2(out_left, out_right) = in_block(i, j);
                }
//            std::cout << "block " << block << " -- finish" <<std::endl;
        }
    }
    
    timer.end();
    
}

// easier code, but slower
template<class Matrix, class SymmGroup>
void reshape2(MultiIndex<SymmGroup> const & midx,
              typename MultiIndex<SymmGroup>::set_id in_set,
              typename MultiIndex<SymmGroup>::set_id out_set,
              block_matrix<Matrix, SymmGroup> const & m1,
              block_matrix<Matrix, SymmGroup> & m2)
{   
    static Timer timer("generic reshape2");
    timer.begin();
    
    
    m2 = block_matrix<Matrix, SymmGroup>();
    
    typedef std::size_t size_t;
    typedef typename SymmGroup::charge charge;
    typedef typename MultiIndex<SymmGroup>::coord_t coord_t;
    
    
    for (int run = 0; run < 2; ++run) {
        bool pretend = (run == 0);
        
        if (!pretend)
            m2.allocate_blocks();
        
        for(index_product_iterator<SymmGroup> it = midx.begin();
            it != midx.end();
            it++)
        {
            coord_t in_left, in_right;
            boost::tie(in_left, in_right) = midx.get_coords(in_set, *it);
            
            if (!m1.has_block(in_left.first, in_right.first))
                continue;
            
            coord_t out_left, out_right;
            boost::tie(out_left, out_right) = midx.get_coords(out_set, *it);
            
            if (pretend)
                m2.reserve_pos(out_left.first, out_right.first,
                               out_left.second, out_right.second);
            else
                m2(out_left, out_right) = m1(in_left, in_right);
        }
    }
    
    timer.end();
}


#endif
