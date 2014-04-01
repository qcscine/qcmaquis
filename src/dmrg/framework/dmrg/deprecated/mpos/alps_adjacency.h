
#ifndef ALPS_ADJACENCY_H
#define ALPS_ADJACENCY_H

#include "adjacency.h"

#include <string>
#include <vector>
#include <algorithm>
#include <map>

namespace adj {
	
	class ALPSAdj : public Adjacency
    {
    public:
		
		ALPSAdj (std::string const &);
		
        std::vector<int> forward(int p) const
		{
			return forward_[p];
		}
        std::vector<int> all(int p) const
		{
			std::vector<int> ret = forward_[p];
			std::copy(backward_[p].begin(), backward_[p].end(), std::back_inserter(ret));
			return ret; 
		}
		
        int size() const
		{
			return size_;
		}
        
        int site_type(int i) { return site_type_map[i]; }
        int bond_type(int i, int j) { return bond_type_map[i][j]; }
		
	private:
		int size_;
		std::vector<std::vector<int> > forward_;
		std::vector<std::vector<int> > backward_;
        
        std::map<int, std::map<int, int> > bond_type_map;
        std::map<int, int> site_type_map;
    };	
	
} // namespace (adj)

#endif