
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
		alps::Parameters parms (ifs);
		graph_type graph (parms);
		
		// storing lattice informations
		size_ = graph.num_sites();
		forward_.resize(size_);
		backward_.resize(size_);
		
		for (graph_type::bond_iterator it=graph.bonds().first; it!=graph.bonds().second; ++it) {
			forward_[graph.vertex_index(graph.source(*it))].push_back(graph.vertex_index(graph.target(*it)));
			backward_[graph.vertex_index(graph.target(*it))].push_back(graph.vertex_index(graph.source(*it)));
		}
	}
}
