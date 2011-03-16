
#ifndef ALPS_ADJACENCY_H
#define ALPS_ADJACENCY_H

#include "adjancency.h"

#include <string>
#include <vector>
#include <algorithm>

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
		
	private:
		int size_;
		std::vector<std::vector<int> > forward_;
		std::vector<std::vector<int> > backward_;		
    };	
	
} // namespace (adj) 

#endif