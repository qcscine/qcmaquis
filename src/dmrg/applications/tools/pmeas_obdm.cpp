/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/block_matrix/detail/alps.hpp"

typedef alps::numeric::matrix<std::complex<double> > Matrix;

#include <alps/hdf5.hpp>

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"

#include "dmrg/models/factory.h"

typedef U1 grp;


template<class Matrix, class SymmGroup>
void measure_correlation_parallel(MPS<Matrix, SymmGroup> const & mps,
                                  const Lattice & lat,
                                  block_matrix<Matrix, SymmGroup> const & identity,
                                  block_matrix<Matrix, SymmGroup> const & fill,
                                  std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops,
                                  std::vector<std::vector<double> >& all_vals,
                                  std::vector<std::vector<std::size_t> >& all_labels)
{
    if (ops.size() != 2)
        throw std::runtime_error("parallel correlation only for 2 operators.");
    
    
    std::vector<double> vals;
    std::vector<std::vector<std::size_t> > num_labels;

#pragma omp parallel for schedule(dynamic)
    for (size_t p1=0; p1<lat.size()-1; ++p1) {
        for (size_t p2=p1+1; p2<lat.size(); ++p2) {
            generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat.size(), identity);
            
            std::vector< std::pair<int, block_matrix<Matrix, SymmGroup> > > op_term;
            op_term.push_back(p1, ops[0].first);
            op_term.push_back(p2, ops[1].first);
            
            generate_mpo::Operator_Term<Matrix, SymmGroup> term;
            term.operators = op_term;
            term.fill_operator = identity;
            mpom.add_term(term);
        
            MPO<Matrix, SymmGroup> mpo = mpom.create_mpo();
            double val = expval(mps, mpo);
            
            std::vector<std::size_t> lab;
            lab.push_bask(p1); lab.push_bask(p2);
            
            vals.push_bask(val);
            num_labels.push_bask(lab);
        }
    }
    int nthread = omp_thread_num();
    all_vals[nthread] = vals;
    all_labels[nthread] = label_strings(lat, num_labels);
}


int main(int argc, char ** argv)
{
    try {
        if (argc != 3)
            throw std::runtime_error("Usage: <parms> <model_parms>");
        
        /// Loading parameters
        std::ifstream param_file(argv[1]);
        if (!param_file)
            throw std::runtime_error("Could not open parameter file.");
        DmrgParameters parms(param_file);
        
        /// Loading model
        std::ifstream model_file(argv[2]);
        if (!model_file)
            throw std::runtime_error("Could not open model file.");
        ModelParameters model_parms(model_file);
        
        /// Parsing model
        boost::shared_ptr<Lattice> lattice;
        boost::shared_ptr<Model<matrix, grp> > model;
        model_parser<matrix, grp>(parms["lattice_library"], parms["model_library"], model_parms,
                                  lattice, model);
        
        Hamiltonian<matrix, grp> H = model->H();
        MPO<matrix, grp> mpo = make_mpo(lattice->size(), H);
        Measurements<matrix, grp> measurements = model->measurements();

        
        /// Initialize & load MPS
        int L = lattice->size();
        MPS<matrix, grp> mps(L);
        {
            alps::hdf5::archive ar(parms["chkpfile"]);
            ar["/state"] >> mps;
        }
        
        
        Measurements<matrix, grp>::mterm_t obdm_meas = measurements.get("Onebody density matrix");
       
        std::vector<std::vector<double> > vals;
        std::vector<std::vector<std::string> > labels;
        int tot_threads;
#pragma omp parallel
        {
#pragma omp single
            {
                tot_threads = omp_num_threads();
                vals.resize(tot_threads);
                labels.resize(tot_threads);
            }
            measure_correlation_parallel(mps, *lattice, measurements.get_identity(), measurements.get_identity(), meas_obdm.operators, vals, labels);
        }
        
        std::vector<double> all_vals;
        std::vector<std::string> all_labels;
        for (int nt=0; nt<tot_threads; ++nt) {
            for (int i=0; i<vals[nt].size(); ++i) {
                all_vals.push_bask( vals[nt][i] );
                all_labels.push_bask( labels[nt][i] );
            }
        }
        
        {
            alps::hdf5::archive ar(parms["resultfile", "w"]);
            ar["/spectrum/results/Onebody density matrix/mean/value"] << all_vals;
            ar["/spectrum/results/Onebody density matrix/labels"] << all_labels;
        }

        
    } catch (std::exception & e) {
        maquis::cerr << "Exception caught:" << std::endl << e.what() << std::endl;
        exit(1);
    }
}
