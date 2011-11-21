
#ifndef MULTI_INDEX_H
#define MULTI_INDEX_H

//#include "dmrg/block_matrix/indexing.h"

#include <vector>
#include <algorithm>
#include <iostream>

#include <boost/operators.hpp>


template <class SymmGroup>
class index_product_iterator
: public boost::forward_iterator_helper<
                                          index_product_iterator<SymmGroup>
                                        , std::vector<std::pair<typename Index<SymmGroup>::charge, std::size_t> >
                                        , std::ptrdiff_t
                                        , std::vector<std::pair<typename Index<SymmGroup>::charge, std::size_t> > *
                                        , std::vector<std::pair<typename Index<SymmGroup>::charge, std::size_t> > &
>
{
public:
    typedef std::size_t index_id;
    typedef Index<SymmGroup> index_t;
    typedef typename index_t::charge charge;    // element of Index
    typedef std::size_t elem_id;                // element of charge
    
    typedef std::pair<charge, elem_id> coord_t; // index inside block_matrix
    typedef std::vector<coord_t> value_type;
    
    typedef std::vector<typename index_t::const_iterator> vec_iterator;
    
    index_product_iterator() : valid(false) { }

    index_product_iterator(std::vector<index_t> const & idx)
    : size(idx.size())
    , begin(size)
    , state(size)
    , end(size)
    , cur_i(size, 0)
    , max_i(size)
    , valid(true)
    {
        for (int i=0; i<size; ++i) {
            begin[i] = idx[i].begin();
            state[i] = idx[i].begin();
            max_i[i] = state[i]->second;
            end[i] = idx[i].end();
        }
    }
    
    value_type operator*() const
    {
        value_type ret(size);
        for (int i=0; i<size; ++i)
            ret[i] = std::make_pair(state[i]->first, cur_i[i]);
        return ret;
    }
    
    void operator++()
    {
        for (int i=size-1; i>=0; --i) {
            cur_i[i]++;
            if (cur_i[i] != max_i[i])
                break;
            
            state[i]++;
            if (state[i] != end[i]) {
                cur_i[i] = 0;
                max_i[i] = state[i]->second;
                break;
            }
            
            if (i == 0) {
                state[i] = end[i];
                valid = false;
                break;
            }
            
            state[i] = begin[i];
            cur_i[i] = 0;
            max_i[i] = state[i]->second;
        }
    }
    
    bool operator==(index_product_iterator<SymmGroup> const & rhs) const
    {
        if (valid != rhs.valid)
            return false;
        if (!valid)
            return true;
        
        return (size == rhs.size) && std::equal(state.begin(), state.end(), rhs.state.begin());
    }
    
    
private:
    
    bool valid;
    std::size_t size;
    
    vec_iterator begin;
    vec_iterator end;
    vec_iterator state;
    
    std::vector<elem_id> cur_i;
    std::vector<std::size_t> max_i;
};

    
    
template <class SymmGroup>
class MultiIndex {
public:
    typedef std::size_t size_t;
    typedef std::size_t index_id;
    typedef Index<SymmGroup> index_t;
    typedef typename index_t::charge charge;    // element of Index
    typedef std::size_t elem_id;                // element of charge
    
    typedef std::pair<charge, elem_id> coord_t; // index inside block_matrix
    
    typedef std::size_t set_id;
    typedef std::vector<coord_t> key_t;
    
    typedef index_product_iterator<SymmGroup> const_iterator;
    
    MultiIndex () { }
    
    index_id insert_index(index_t const & idx)
    {
        idx_.push_back(idx);
        return idx_.size()-1;
    }
    
    index_t const& index(index_id i) const
    { return idx_[i]; }
    
    size_t size() const
    { return idx_.size(); }
    
    set_id create_set(std::vector<std::pair<index_id, bool> > const & vec_left,
                      std::vector<std::pair<index_id, bool> > const & vec_right);
    
    void remove_set(set_id s)
    {
        left_keys.erase(left_keys.begin() + s);
        left_vals.erase(left_vals.begin() + s);
        set_left.erase(set_left.begin() + s);
        right_keys.erase(right_keys.begin() + s);
        right_vals.erase(right_vals.begin() + s);
        set_right.erase(set_right.begin() + s);
    }
    
    void clear_sets()
    {
        left_keys.clear();
        left_vals.clear();
        left_sizes.clear();
        set_left.clear();
        right_keys.clear();
        right_vals.clear();
        right_sizes.clear();
        set_right.clear();
    }
    
    void clear()
    {
        clear_sets();
        idx_.clear();
    }
    
    
    size_t left_size(set_id s, charge const & c) const
    {
        assert( left_sizes[s].count(c) > 0 );
        return left_sizes[s].find(c)->second;
    }
    size_t right_size(set_id s, charge const & c) const
    {
        assert( right_sizes[s].count(c) > 0 );
        return right_sizes[s].find(c)->second;
    }
    
    // TODO: cache the last search key to avoid find()
    
    // key --> coord
    inline coord_t const& get_left_coord(set_id s, key_t const& key) const
    {
        return key_to_val(key, left_keys[s], left_vals[s]);
    }
    inline coord_t const& get_right_coord(set_id s, key_t const& key) const
    {
        return key_to_val(key, right_keys[s], right_vals[s]);
    }
    std::pair<coord_t, coord_t> get_coords(set_id s, key_t const& key) const
    {
        key_t left_k, right_k;
        for (int i=0; i<set_left[s].size(); ++i)
            left_k.push_back(key[ set_left[s][i].first ]);
        for (int i=0; i<set_right[s].size(); ++i)
            right_k.push_back(key[ set_right[s][i].first ]);
                
        return std::make_pair( get_left_coord(s, left_k), get_right_coord(s, right_k) );
    }
    
    // coord --> key
    inline key_t const& get_left_key(set_id s, coord_t const& ci) const
    {
        return val_to_key(ci, left_keys[s], left_vals[s]);
    }
    inline key_t const& get_right_key(set_id s, coord_t const& ci) const
    {
        return val_to_key(ci, right_keys[s], right_vals[s]);
    }
    
    // coord --> coord
    std::pair<coord_t, coord_t> convert_coords(set_id set1,
                                               coord_t const & coord1_left,
                                               coord_t const & coord1_right,
                                               set_id set2) const
    {
        key_t key1_left = get_left_key(set1, coord1_left);
//        std::cout << "got left key" << std::endl;
        key_t key1_right = get_right_key(set1, coord1_right);
//        std::cout << "got right key" << std::endl;

        assert( key1_left.size() == set_left[set1].size() );
        assert( key1_right.size() == set_right[set1].size() );

        key_t key(idx_.size());
        for (int i=0; i<set_left[set1].size(); ++i) {
//            std::cout << "key index = " << set_left[set1][i].first << std::endl;
            key[ set_left[set1][i].first ] = key1_left[i];
        }
//        std::cout << "left key --> key" << std::endl;
        for (int i=0; i<set_right[set1].size(); ++i)
            key[ set_right[set1][i].first ] = key1_right[i];
//        std::cout << "right key --> key" << std::endl;

        return get_coords(set2, key);
        
    }
    inline std::pair<coord_t, coord_t> convert_coords(set_id set1, std::pair<coord_t, coord_t> const & coords1,
                                               set_id set2) const
    { return convert_coords(set1, coords1.first, coords1.second, set2); }
    
    
    const_iterator begin() const
    { return const_iterator(idx_); }
    
    const_iterator end() const
    { return const_iterator(); }
    
    
private:
    coord_t const& key_to_val(key_t const& key, std::vector<key_t> const& keys,
                              std::vector<coord_t> const& vals) const
    {
        assert( std::count(keys.begin(), keys.end(), key) > 0 );
        size_t pos = std::find(keys.begin(), keys.end(), key)-keys.begin();
        return vals[pos];
    }
    
    key_t const& val_to_key(coord_t const& ci, std::vector<key_t> const& keys,
                            std::vector<coord_t> const& vals) const
    {
        assert( std::count(vals.begin(), vals.end(), ci) > 0 );
        size_t pos = std::find(vals.begin(), vals.end(), ci)-vals.begin();
        return keys[pos];
    }
    
    
    
    std::vector<index_t> idx_;
    
    std::vector<std::vector<std::pair<index_id, bool> > > set_left, set_right;
    std::vector<std::vector<key_t> > left_keys, right_keys;
    std::vector<std::vector<coord_t> > left_vals, right_vals;
    std::vector<std::map<charge, size_t> > left_sizes, right_sizes;
    
};


template <class SymmGroup>
typename MultiIndex<SymmGroup>::set_id
MultiIndex<SymmGroup>::create_set(std::vector<std::pair<index_id, bool> > const & vec_left,
                                  std::vector<std::pair<index_id, bool> > const & vec_right)
{
    {
        std::vector<index_t> b;
        for(int i=0; i<vec_left.size(); ++i)
            b.push_back( idx_[vec_left[i].first] );
        
        std::map<charge, size_t> block_begins;
        std::vector<key_t> keys_;
        std::vector<coord_t> vals_;
        
        for (const_iterator it = const_iterator(b);
             it != const_iterator();
             it++)
        {
            charge c = SymmGroup::IdentityCharge;
            for (int i=0; i<b.size(); ++i)
                c = (vec_left[i].second) ? SymmGroup::fuse(c, (*it)[i].first) : SymmGroup::fuse(c, -(*it)[i].first);
            keys_.push_back(*it);
            vals_.push_back(std::make_pair(c, block_begins[c]));
            block_begins[c]++;
        }
        
        left_keys.push_back(keys_);
        left_vals.push_back(vals_);
        left_sizes.push_back(block_begins);
    }
    
    {
        std::vector<index_t> b;
        for(int i=0; i<vec_right.size(); ++i)
            b.push_back( idx_[vec_right[i].first] );
        
        std::map<charge, size_t> block_begins;
        std::vector<key_t> keys_;
        std::vector<coord_t> vals_;
        
        for (const_iterator it = const_iterator(b);
             it != const_iterator();
             it++)
        {
            charge c = SymmGroup::IdentityCharge;
            for (int i=0; i<b.size(); ++i)
                c = (vec_right[i].second) ? SymmGroup::fuse(c, (*it)[i].first) : SymmGroup::fuse(c, -(*it)[i].first);
            keys_.push_back(*it);
            vals_.push_back(std::make_pair(c, block_begins[c]));
            block_begins[c]++;
        }
        
        right_keys.push_back(keys_);
        right_vals.push_back(vals_);
        right_sizes.push_back(block_begins);
    }
    
    set_left.push_back(vec_left);
    set_right.push_back(vec_right);
    return left_keys.size() - 1;
}

// ostreams
//template <class SymmGroup>
//std::ostream& operator<< (std::ostream& os, std::pair<typename SymmGroup::charge, std::size_t> const& p)
//{
//    os << "(" << p.first << " : " << p.second << ")";
//    return os;
//}
//
//template <class SymmGroup>
//std::ostream& operator<< (std::ostream& os, typename index_product_iterator<SymmGroup>::value_type const& v)
//{
//    //std::copy(v.begin(), v.end(), std::ostream_iterator<std::pair<symm::charge, std::size_t> >(os, " "));
//    for (int i=0; i<v.size(); ++i)
//        os << v[i] << " ";
//    return os;
//}
//
//template <class SymmGroup>
//std::ostream& operator<< (std::ostream& os, std::pair<typename MultiIndex<SymmGroup>::coord_t, typename MultiIndex<SymmGroup>::coord_t> const& p)
//{
//    os << p.first << ", " << p.second;
//    return os;
//}



#endif
