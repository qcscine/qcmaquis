#ifndef ALPS_LATTICE_H
#define ALPS_LATTICE_H

#include "lattice.h"

#include <vector>
#include <map>
#include <algorithm>
#include <sstream>

#include <alps/parameter.h>
#include <alps/lattice.h>


class ALPSLattice : public Lattice
{

public:
    typedef Lattice::pos_t pos_t;
    typedef alps::graph_helper<> graph_type;
    typedef graph_type::site_descriptor site_descriptor;
    typedef graph_type::site_iterator site_iterator;
    
    ALPSLattice (std::ifstream & ifs) : 
    parms(ifs),
    graph(parms)
    {
        // storing lattice informations
        forward_.resize(size());
        backward_.resize(size());
                
        for (graph_type::bond_iterator it=graph.bonds().first; it!=graph.bonds().second; ++it) {
            graph_type::size_type s, t;
            s = graph.vertex_index(graph.source(*it));
            t = graph.vertex_index(graph.target(*it));
            
            forward_[s].push_back(t);
            backward_[t].push_back(s);
            
            bond_index_map[s][t] = graph.edge_index(*it);
            bond_index_map[t][s] = graph.edge_index(*it);
        }
    }
    
    std::vector<pos_t> forward(pos_t p) const
    {
        return forward_[p];
    }
    std::vector<pos_t> all(pos_t p) const
    {
        std::vector<pos_t> ret = forward_[p];
        std::copy(backward_[p].begin(), backward_[p].end(), std::back_inserter(ret));
        return ret; 
    }
    
    pos_t size() const
    {
        return graph.num_sites();
    }
    
    boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const
    {
        if (property == "label" && pos.size() == 1)
            return boost::any( alps::site_label(graph.graph(), graph.site(pos[0])) );
        else if (property == "label" && pos.size() == 2)
            return boost::any( alps::bond_label(graph.graph(), graph.bond(bond_index_map[pos[0]][pos[1]])) );
        else if (property == "type" && pos.size() == 1)
            return boost::any( static_cast<int>(graph.site_type(graph.site(pos[0]))) );
        else if (property == "type" && pos.size() == 2)
            return boost::any( static_cast<int>(graph.bond_type(graph.bond(bond_index_map[pos[0]][pos[1]]))) );
        else if (property == "wraps_pbc" && pos.size() == 2)
            return boost::any( static_cast<bool>(boost::get(alps::boundary_crossing_t(),
                                                            graph.graph(),
                                                            graph.bond(bond_index_map[pos[0]][pos[1]]))) );
        else {
            std::ostringstream ss;
            ss << "No property '" << property << "' with " << pos.size() << " points implemented."; 
            throw std::runtime_error(ss.str());
            return boost::any();
        }
    }

    
private:
    alps::Parameters parms;
    graph_type graph;
    std::vector<std::vector<pos_t> > forward_;
    std::vector<std::vector<pos_t> > backward_;
    // TODO: find better solution instead of this (operator[] is non-const)
    mutable std::map<pos_t, std::map<pos_t, pos_t> > bond_index_map;
};


#endif