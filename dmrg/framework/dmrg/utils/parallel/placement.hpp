/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef PLACEMENT_HPP
#define PLACEMENT_HPP

namespace parallel {

inline void make_consistent(std::vector<int>& p_new, int b2, std::vector<int>& e_new,
                            const std::vector<int>& p_old, int b1, std::vector<int>& e_old)
{
    auto it_old = std::find(e_old.begin(), e_old.end(), b1);
    auto it_new = std::find(e_new.begin(), e_new.end(), b2);
    if (it_old != e_old.end() && it_new != e_new.end()) {
        auto it_inner = std::find(p_new.begin(), p_new.end(), p_old[b1]);
        if (it_inner != p_new.end()) {
            e_old.push_back(b1);
            std::replace(p_new.begin(), p_new.end(), p_old[b1], -1);
        }
        else if(p_new[b2] != -1) {
            e_new.push_back(b2);
            p_new[b2] = -1;
        } else {
            p_new[b2] = p_old[b1];
        }
    }
}

inline void make_complete(std::vector<int>& placement) {
    int iToAdd = 0;
    for(int i = 0; i < placement.size(); i++) {
        if(placement[i] == -1) {
            bool added = false;
            while (!added) {
                auto itFind = std::find(placement.begin(), placement.end(), iToAdd);
                if(itFind == placement.end()) {
                    placement[i] = iToAdd;
                    added = true;
                }
                else {
                    iToAdd += 1;
                }
            }
        }
    }
}

template<class Matrix, class SymmGroup>
std::vector<int> get_left_placement(const MPOTensor<Matrix, SymmGroup>& mpo, const std::vector<int>& placement_l_old, const std::vector<int>& placement_r) {
    std::vector<int> ts_exceptions_l, ts_exceptions_r;
    std::vector<int> placement_l(placement_l_old.size(), -1);
    ts_exceptions_l.reserve(placement_l_old.size());
    ts_exceptions_r.reserve(placement_r.size());
    for(int b1 = 0; b1 < placement_l_old.size(); b1++)
        for(int b2 = 0; b2 < placement_r.size(); b2++) {
            if(mpo.has(b1,b2)) {
                make_consistent(placement_l, b1, ts_exceptions_l,
                                placement_r, b2, ts_exceptions_r);
            }
    }
    make_complete(placement_l);
    mpo.exceptions_l = ts_exceptions_l;
    mpo.exceptions_r = ts_exceptions_r;
    return placement_l;
}

template<class Matrix, class SymmGroup>
std::vector<int> get_right_placement(const MPOTensor<Matrix, SymmGroup>& mpo, const std::vector<int>& placement_l, const std::vector<int>& placement_r_old) {
    std::vector<int> ts_exceptions_l, ts_exceptions_r;
    std::vector<int> placement_r(placement_r_old.size(), -1);
    ts_exceptions_l.reserve(placement_l.size());
    ts_exceptions_r.reserve(placement_r_old.size());
    for(int b1 = 0; b1 < placement_l.size(); b1++)
    for(int b2 = 0; b2 < placement_r_old.size(); b2++){
        if(mpo.has(b1,b2)){
            make_consistent(placement_r, b2, ts_exceptions_r, 
                            placement_l, b1, ts_exceptions_l);
        }
    }
    make_complete(placement_r);
    mpo.exceptions_l = ts_exceptions_l;
    mpo.exceptions_r = ts_exceptions_r;
    return placement_r;
}

template<class Matrix, class SymmGroup>
std::vector<std::vector<int> > construct_placements(const MPO<Matrix, SymmGroup>& mpo) {
    std::vector<std::vector<int> > placements(mpo.length()+1);
    std::vector<std::pair<std::vector<int>, std::vector<int> > > exceptions(mpo.length()+1);
    placements[0].push_back(0); // left_[0] has only 1 element
    for(int s = 0; s < mpo.length(); s++){
        placements[s+1].resize(mpo[s].col_dim(), -1);

        for(int b1 = 0; b1 < placements[s].size(); b1++)
        for(int b2 = 0; b2 < placements[s+1].size(); b2++){
            if(mpo[s].has(b1,b2)){
                make_consistent(placements[s+1], b2, exceptions[s+1].second, placements[s], b1, exceptions[s].first);
            }
        }
        make_complete(placements[s+1]);
        mpo[s].placement_l = placements[s];
        mpo[s].placement_r = placements[s+1];
        mpo[s].exceptions_l = exceptions[s].first;
        mpo[s].exceptions_r = exceptions[s+1].second;
    }
    return placements;
}

}

#endif


