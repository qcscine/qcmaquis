/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "mp_tensors/mpo.h"
#include "mp_tensors/mps_mpo_ops.h"

#include <sstream>

#include <alps/hdf5.hpp>

namespace meas_detail {
    
	inline std::vector<std::string> label_strings (const Lattice& lat, const std::vector<std::vector<std::size_t> >& labels)
	{
		std::vector<std::string> ret(labels.size());
		std::vector<std::string>::iterator ot = ret.begin();
		for (std::vector<std::vector<std::size_t> >::const_iterator it = labels.begin();
			 it != labels.end(); ++it)
		{
			std::ostringstream oss;
			for (std::vector<std::size_t>::const_iterator it2 = it->begin(); it2 != it->end(); ++it2) {
				oss << lat.get_prop<std::string>("label", *it2);
				if (it2 + 1 != it->end())
					oss << " -- ";
			}
			*(ot++) = oss.str();
		}
		assert( ot == ret.end() );
		return ret;
	}
    
    
	template<class Matrix, class SymmGroup>
	void measure_local(MPS<Matrix, SymmGroup> & mps,
                       const Lattice & lat,
                       block_matrix<Matrix, SymmGroup> const & identity,
                       block_matrix<Matrix, SymmGroup> const & op,
                       alps::hdf5::oarchive & ar,
                       std::string base_path)
	{
		std::vector<double> vals;
        std::vector<std::string> labels;
		for (std::size_t p = 0; p < lat.size(); ++p)
        {
            generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat.size(), identity);
            generate_mpo::Operator_Term<Matrix, SymmGroup> term;
            term.operators.push_back( std::make_pair(p, op) );
            mpom.add_term(term);
            MPO<Matrix, SymmGroup> mpo = mpom.create_mpo();
            
            double val = expval(mps, mpo, 0);
            vals.push_back(val);
            labels.push_back(lat.get_prop<std::string>("label", p));
        }
        
	    ar << alps::make_pvp(base_path + std::string("/mean/value"), vals);
	    ar << alps::make_pvp(base_path + std::string("/labels"), labels);
	}
    
	template<class Matrix, class SymmGroup>
	void measure_average(MPS<Matrix, SymmGroup> & mps,
                         const Lattice & lat,
                         block_matrix<Matrix, SymmGroup> const & identity,
                         block_matrix<Matrix, SymmGroup> const & op,
                         alps::hdf5::oarchive & ar,
                         std::string base_path)
	{
        generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat.size(), identity);
		for (std::size_t p = 0; p < lat.size(); ++p)
        {
            generate_mpo::Operator_Term<Matrix, SymmGroup> term;
            term.operators.push_back( std::make_pair(p, op) );
            mpom.add_term(term);
        }
        MPO<Matrix, SymmGroup> mpo = mpom.create_mpo();
        double val = expval(mps, mpo, 0);
        
	    ar << alps::make_pvp(base_path + std::string("/mean/value"), val);
	}
    
	template<class MAKER, class Matrix, class SymmGroup>
	void measure_correlation(MPS<Matrix, SymmGroup> & mps,
							 const Lattice & lat,
							 block_matrix<Matrix, SymmGroup> const & identity,
							 block_matrix<Matrix, SymmGroup> const & fill,
							 std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops,
							 alps::hdf5::oarchive & ar,
							 std::string base_path,
							 bool half=false)
	{
	    std::vector<double> dc;
	    std::vector<std::string> labels;
	    if (half) {
	    	for (size_t p = 0; p < lat.size()-(ops.size()-1); ++p) {
	    		MAKER dcorr(mps.length(), identity, fill,
                            ops, p);
	    		MPO<Matrix, SymmGroup> mpo = dcorr.create_mpo();
                
	    		std::cout << dcorr.description() << std::endl;
                
	    		std::vector<double> dct = multi_expval(mps, mpo);
	    		std::copy(dct.begin(), dct.end(), std::back_inserter(dc));
                
	    		std::vector<std::string> lbt = label_strings(lat, dcorr.numeric_labels());
	    		std::copy(lbt.begin(), lbt.end(), std::back_inserter(labels));
	    	}
	    } else {
	    	MAKER dcorr(mps.length(), identity, fill,
                        ops);
	        MPO<Matrix, SymmGroup> mpo = dcorr.create_mpo();
            
    		std::cout << dcorr.description() << std::endl;
            
	        dc = multi_expval(mps, mpo);
	        labels = label_strings(lat, dcorr.numeric_labels());
	    }
	    ar << alps::make_pvp(base_path + std::string("/mean/value"), dc);
	    ar << alps::make_pvp(base_path + std::string("/labels"), labels);
	}
    
    
} // namespace
