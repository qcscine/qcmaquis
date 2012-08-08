#include <iostream>

#include <vector>
#include <iterator>
#include <algorithm>

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/multi_index.h"
#include "dmrg/block_matrix/detail/alps.hpp"
typedef alps::numeric::matrix<double> Matrix;

typedef U1 symm;

std::ostream& operator<< (std::ostream& os, std::pair<symm::charge, std::size_t> const& p)
{
	os << "(" << p.first << " : " << p.second << ")";
	return os;
}

std::ostream& operator<< (std::ostream& os, index_product_iterator<symm>::value_type const& v)
{
	//std::copy(v.begin(), v.end(), std::ostream_iterator<std::pair<symm::charge, std::size_t> >(os, " "));
	for (int i=0; i<v.size(); ++i)
		os << v[i] << " ";
	return os;
}

std::ostream& operator<< (std::ostream& os, std::pair<MultiIndex<symm>::coord_t, MultiIndex<symm>::coord_t> const& p)
{
	os << p.first << ", " << p.second;
	return os;
}

typedef std::vector<std::pair<symm::charge, std::size_t> > access_t;
typedef std::vector<std::pair<std::pair<symm::charge, std::size_t>, std::pair<symm::charge, std::size_t> > > map_val_t;
typedef std::map<access_t,  map_val_t> ret_map;

ret_map test_multi_index(Index<symm> const & ix1, Index<symm> const & ix2, Index<symm> const & phys) {
	typedef std::vector< std::pair<MultiIndex<symm>::index_id, bool> > set_descriptor_t;
	
	// construct and fill idx_
	MultiIndex<symm> midx;
	MultiIndex<symm>::index_id alpha = midx.insert_index(ix1);
	MultiIndex<symm>::index_id beta = midx.insert_index(ix2);
	MultiIndex<symm>::index_id s1 = midx.insert_index(phys);
	MultiIndex<symm>::index_id s2 = midx.insert_index(phys);

	// create sets
	set_descriptor_t set1_left(3), set1_right(1);
	set1_left[0] = std::make_pair(s1, true); set1_left[1] = std::make_pair(s2, true); set1_left[2] = std::make_pair(alpha, true);
	set1_right[0] = std::make_pair(beta, true);
	MultiIndex<symm>::set_id set1 = midx.create_set(set1_left, set1_right);
	
	set_descriptor_t set2_left(1), set2_right(3);
	set2_left[0] = std::make_pair(alpha, true);
	set2_right[0] = std::make_pair(s1, false); set2_right[1] = std::make_pair(s2, false); set2_right[2] = std::make_pair(beta, true); 
	MultiIndex<symm>::set_id set2 = midx.create_set(set2_left, set2_right);

    set_descriptor_t set3_left(2), set3_right(2);
	set3_left[0] = std::make_pair(s1, true); set3_left[1] = std::make_pair(alpha, true);
	set3_right[0] = std::make_pair(s2, false); set3_right[1] = std::make_pair(beta, true);
	MultiIndex<symm>::set_id set3 = midx.create_set(set3_left, set3_right);


    ret_map ret;
	for(index_product_iterator<symm> it = midx.begin();
		it != midx.end();
		it++)
	{
		maquis::cout << *it << " = " << midx.get_coords(set1, *it)     // left pairing
                            << " or " << midx.get_coords(set2, *it)    // right pairing
                            << " or " << midx.get_coords(set3, *it)    // both pairing
                            << std::endl;
        ret[*it].push_back(midx.get_coords(set1, *it));
        ret[*it].push_back(midx.get_coords(set2, *it));
        ret[*it].push_back(midx.get_coords(set3, *it));
	}

    return ret;
}


ret_map test_usual_l(Index<symm> const & left_i, Index<symm> const & right_i, Index<symm> const & physical_i) {
    typedef std::size_t size_t;
    typedef symm::charge charge;
    
    Index<symm> phys2_i = physical_i*physical_i;
    ProductBasis<symm> phys_pb(physical_i, physical_i);
    ProductBasis<symm> left_both(physical_i, left_i);
    ProductBasis<symm> right_both(physical_i, right_i,
                                  boost::lambda::bind(static_cast<charge(*)(charge, charge)>(symm::fuse),
                                                      -boost::lambda::_1, boost::lambda::_2));
    ProductBasis<symm> left_p(phys2_i, left_i);
    ProductBasis<symm> right_p(phys2_i, right_i,
                               boost::lambda::bind(static_cast<charge(*)(charge, charge)>(symm::fuse),
                                                   -boost::lambda::_1, boost::lambda::_2));
    
    ret_map ret;
    for (size_t l = 0; l < left_i.size(); ++l)
        for (size_t r = 0; r < right_i.size(); ++r)
            for (size_t s1 = 0; s1 < physical_i.size(); ++s1)
                for (size_t s2 = 0; s2 < physical_i.size(); ++s2)
                {
                    charge s_charge = symm::fuse(physical_i[s1].first, physical_i[s2].first);
                    size_t s_out = phys2_i.find(s_charge);
                    
                    charge bp_l_charge = symm::fuse(physical_i[s1].first, left_i[l].first);    // both paired,  left
                    charge bp_r_charge = symm::fuse(-physical_i[s2].first, right_i[r].first);  // both paired,  right
                    charge lp_l_charge = symm::fuse(s_charge, left_i[l].first);                // left paired,  left
                    charge lp_r_charge = right_i[r].first;                                     // left paired,  right
                    charge rp_l_charge = left_i[l].first;                                      // right paired, left
                    charge rp_r_charge = symm::fuse(-s_charge, right_i[r].first);              // right paired, right
                    
                    size_t bp_left_offset = left_both(physical_i[s1].first, left_i[l].first);
                    size_t bp_right_offset = right_both(physical_i[s2].first, right_i[r].first);
                    size_t lp_left_offset = left_p(s_charge, left_i[l].first);
                    size_t lp_right_offset = 0;
                    size_t rp_left_offset = 0;
                    size_t rp_right_offset = right_p(s_charge, right_i[r].first);
                    size_t out_phys_offset = phys_pb(physical_i[s1].first, physical_i[s2].first);
                                        
                    for (size_t ss1 = 0; ss1 < physical_i[s1].second; ++ss1)
                        for (size_t ss2 = 0; ss2 < physical_i[s2].second; ++ss2)
                        {
                            size_t ss_out = out_phys_offset + ss1*physical_i[s2].second + ss2;
                            for (size_t rr = 0; rr < right_i[r].second; ++rr)
                                for (size_t ll = 0; ll < left_i[l].second; ++ll)
                                {
                                    access_t aa(4);
                                    aa[0] = std::make_pair(left_i[l].first, ll);
                                    aa[1] = std::make_pair(right_i[r].first, rr);
                                    aa[2] = std::make_pair(physical_i[s1].first, ss1);
                                    aa[3] = std::make_pair(physical_i[s2].first, ss2);
                                    
                                    map_val_t vv(3);
                                    // left pairing
                                    vv[0] = std::make_pair(std::make_pair(lp_l_charge, lp_left_offset + ss_out*left_i[l].second + ll),
                                                           std::make_pair(lp_r_charge, lp_right_offset + rr));
                                    // right pairing
                                    vv[1] = std::make_pair(std::make_pair(rp_l_charge, rp_left_offset + ll),
                                                           std::make_pair(rp_r_charge, rp_right_offset + ss_out*right_i[r].second + rr));
                                    // both pairing
                                    vv[2] = std::make_pair(std::make_pair(bp_l_charge, bp_left_offset + ss1*left_i[l].second+ll),
                                                           std::make_pair(bp_r_charge, bp_right_offset + ss2*right_i[r].second+rr));
                                    
                                    maquis::cout << aa << " = " << vv[0]
                                                       << " or " << vv[1]
                                                       << " or " << vv[2]
                                                       << std::endl;
                                    
                                    ret[aa] = vv;
                                }
                        }
                }
    return ret;
}

ret_map test_usual_r(Index<symm> const & left_i, Index<symm> const & right_i, Index<symm> const & physical_i) {
    typedef std::size_t size_t;
    typedef symm::charge charge;
    
    Index<symm> pyleft_i = physical_i*left_i;
    Index<symm> pyright_i = adjoin(physical_i)*right_i;
    ProductBasis<symm> left_both(physical_i, left_i);
    ProductBasis<symm> right_both(physical_i, right_i,
                                  boost::lambda::bind(static_cast<charge(*)(charge, charge)>(symm::fuse),
                                                      -boost::lambda::_1, boost::lambda::_2));
    ProductBasis<symm> left_p(physical_i, pyleft_i);
    ProductBasis<symm> right_p(physical_i, pyright_i,
                               boost::lambda::bind(static_cast<charge(*)(charge, charge)>(symm::fuse),
                                                   -boost::lambda::_1, boost::lambda::_2));
    
    ret_map ret;
    for (size_t l = 0; l < left_i.size(); ++l)
        for (size_t r = 0; r < right_i.size(); ++r)
            for (size_t s1 = 0; s1 < physical_i.size(); ++s1)
                for (size_t s2 = 0; s2 < physical_i.size(); ++s2)
                {
                    charge pyl_charge = symm::fuse(physical_i[s2].first, left_i[l].first);
                    size_t pyl = pyleft_i.find(pyl_charge);
                    
                    charge pyr_charge = symm::fuse(-physical_i[s2].first, right_i[r].first);
                    size_t pyr = pyright_i.find(pyr_charge);
                    
                    
                    charge bp_l_charge = symm::fuse(physical_i[s1].first, left_i[l].first);    // both paired,  left
                    charge bp_r_charge = symm::fuse(-physical_i[s2].first, right_i[r].first);  // both paired,  right
                    charge lp_l_charge = symm::fuse(physical_i[s1].first, pyl_charge);                // left paired,  left
                    charge lp_r_charge = right_i[r].first;                                     // left paired,  right
                    charge rp_l_charge = left_i[l].first;                                      // right paired, left
                    charge rp_r_charge = symm::fuse(-physical_i[s1].first, pyr_charge);              // right paired, right
                    
                    size_t bp_left_offset = left_both(physical_i[s1].first, left_i[l].first);
                    size_t bp_right_offset = right_both(physical_i[s2].first, right_i[r].first);
                    size_t lp_left_offset = left_p(physical_i[s1].first, pyl_charge);
                    size_t lp_right_offset = 0;
                    size_t rp_left_offset = 0;
                    size_t rp_right_offset = right_p(physical_i[s1].first, pyr_charge);
                    size_t pyl_offset = left_both(physical_i[s2].first, left_i[l].first);
                    size_t pyr_offset = right_both(physical_i[s2].first, right_i[r].first);
                    
                    for (size_t ss1 = 0; ss1 < physical_i[s1].second; ++ss1)
                        for (size_t ss2 = 0; ss2 < physical_i[s2].second; ++ss2)
                        {
                            for (size_t rr = 0; rr < right_i[r].second; ++rr)
                                for (size_t ll = 0; ll < left_i[l].second; ++ll)
                                {
                                    size_t pyll = pyl_offset + ss2*left_i[l].second + ll;
                                    size_t pyrr = pyr_offset + ss2*right_i[r].second + rr;
                                    
                                    access_t aa(4);
                                    aa[0] = std::make_pair(left_i[l].first, ll);
                                    aa[1] = std::make_pair(right_i[r].first, rr);
                                    aa[2] = std::make_pair(physical_i[s1].first, ss1);
                                    aa[3] = std::make_pair(physical_i[s2].first, ss2);
                                    
                                    map_val_t vv(3);
                                    // left pairing
                                    vv[0] = std::make_pair(std::make_pair(lp_l_charge, lp_left_offset + ss1*pyleft_i[pyl].second + pyll),
                                                           std::make_pair(lp_r_charge, lp_right_offset + rr));
                                    // right pairing
                                    vv[1] = std::make_pair(std::make_pair(rp_l_charge, rp_left_offset + ll),
                                                           std::make_pair(rp_r_charge, rp_right_offset + ss1*pyright_i[pyr].second + pyrr));
                                    // both pairing
                                    vv[2] = std::make_pair(std::make_pair(bp_l_charge, bp_left_offset + ss1*left_i[l].second+ll),
                                                           std::make_pair(bp_r_charge, bp_right_offset + ss2*right_i[r].second+rr));
                                    
                                    maquis::cout << aa << " = " << vv[0]
                                    << " or " << vv[1]
                                    << " or " << vv[2]
                                    << std::endl;
                                    
                                    ret[aa] = vv;
                                }
                        }
                }
    return ret;
}


int main () {
	
	Index<symm> ix1;
	ix1.insert(std::make_pair(0, 1));
	ix1.insert(std::make_pair(1, 1));

	Index<symm> ix2;
	ix2.insert(std::make_pair(0, 2));
	ix2.insert(std::make_pair(1, 1));
	
	Index<symm> phys;
	phys.insert(std::make_pair(0, 1));
	phys.insert(std::make_pair(1, 2));
	phys.insert(std::make_pair(2, 1));
	

	maquis::cout << "*** MultiIndex" << std::endl;
	ret_map multi_res = test_multi_index(ix1, ix2, phys);

	maquis::cout << "*** Usual reshape: iix1 = phys * (phys * ix1)" << std::endl;
	ret_map right_res = test_usual_r(ix1, ix2, phys);

	maquis::cout << "*** Usual reshape: iix1 = (phys * phys) * ix1" << std::endl;
	ret_map left_res = test_usual_l(ix1, ix2, phys);
    
    
    bool multi_right=true, multi_left=true, right_left=true;
    
    for (ret_map::const_iterator it=multi_res.begin(); it!= multi_res.end(); ++it)
    {
        if (multi_right)
            multi_right = std::equal(it->second.begin(), it->second.end(), right_res[it->first].begin());
        if (multi_left)
            multi_left = std::equal(it->second.begin(), it->second.end(), left_res[it->first].begin());
        if (right_left)
            right_left = std::equal(right_res[it->first].begin(), right_res[it->first].end(), left_res[it->first].begin());
    }
    
    maquis::cout << "multi_index is equal to right_first? " << ( (multi_right) ? "yes" : "no" ) << std::endl;
    maquis::cout << "multi_index is equal to left_first? " << ( (multi_left) ? "yes" : "no" ) << std::endl;
    maquis::cout << "right_first is equal to left_first? " << ( (right_left) ? "yes" : "no" ) << std::endl;
}

