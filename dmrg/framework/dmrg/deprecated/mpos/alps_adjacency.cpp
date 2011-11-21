
#include "alps_adjacency.h"

#include <fstream>

#include <alps/parameter.h>
#include <alps/lattice.h>

namespace adj {

	ALPSAdj::ALPSAdj (std::string const & parms_file)
	{
		typedef alps::graph_helper<> graph_type;
		
		// initializing the lattice
		std::ifstream ifs(parms_file.c_str());
        if (!ifs)
            throw std::runtime_error("ALPS lattice file not found.");
        
		alps::Parameters parms (ifs);
		graph_type graph (parms);
		
		// storing lattice informations
		size_ = graph.num_sites();
		forward_.resize(size_);
		backward_.resize(size_);
		
		for (graph_type::bond_iterator it=graph.bonds().first; it!=graph.bonds().second; ++it) {
            graph_type::size_type s, t;
            s = graph.vertex_index(graph.source(*it));
            t = graph.vertex_index(graph.target(*it));
            
			forward_[s].push_back(t);
			backward_[t].push_back(s);
            
            bond_type_map[s][t] = graph.bond_type(*it);
            bond_type_map[t][s] = graph.bond_type(*it);
		}
        
        for (graph_type::site_iterator it = graph.sites().first; it != graph.sites().second; ++it)
            site_type_map[graph.vertex_index(*it)] = graph.site_type(*it);
	}
}
