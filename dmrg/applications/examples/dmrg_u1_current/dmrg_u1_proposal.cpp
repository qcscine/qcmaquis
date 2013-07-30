/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *                            Bela Bauer <bauerb@comp-phys.org>
 *
 *****************************************************************************/

#include "maquis/dmrg.hpp"
#include <alps/numeric/matrix.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>

typedef alps::numeric::matrix<double> matrix;
typedef U1 symm;


int main(int argc, char ** argv)
{
    try {
      // use exceptions unstead of exit(1), same below
      if (argc != 3)
        throw std::runtime_error("Usage: <parms> <model_parms>");
      
      maquis::parameters parms(argv[1],maquis::dmrg::default_parms);  
      maquis::parameters model_parms(argv[2],maquis::model::default_parms);
      maquis::storage::setup(parms);

      /// Parsing model
      maquis::lattice lattice(model_parms);
      maquis::model<matrix, symm> model(model_parms,lattice);
      
      maquis::dmrg::mpo<matrix, symm> mpo = model.mpo();
      
      /// Initialize MPS
      maquis::dmrg::mps<matrix, symm> mps = model.intialize_mps(parms["initial_bond_dimension"],parms);
                                 //model.physical_dimension(), model.total_charge(),
      
      /// Initialize optimizer /* , ... restart(sweep, -starting site) */)
      maquis::dmrg::optimizer<matrix, symm, maquis::single_site, maquis::storage::disk> optimizer(mps, mpo, parms);
      // could be implemented as template <class M, class S, class P, class D> class optimizer : P::apply<M,S,D>::type;
      
      /// sweep
      int sweeps = parms["nsweeps"];
      for (int sweep=0; sweep<sweeps; ++sweep) {
        // use exceptions if time limit is exceeded. time limit was in parms
        optimizer.sweep();
        // alternative: maquis::dmrg::single_site_sweep(mps,mpo,parms);
        storage::archive ar(parms["resultfile"].str(), "w");
        std::string iteration_path = std::string("/simulation/sweeps/")+sweep; // maybe boost::lexical_cast
        ar[iteration_path + "/parameters"] << parms;
        ar[iteration_path + "/parameters"] << model_parms;
        ar[iteration_path + "/results"] << optimizer.iteration_results();
      }
      // drop it since the optimizer should reference the mps
      //  mps = optimizer.get_current_mps();
      
      /// Measurement on final MPS
      maquis::measurements results = maquis::measure(mps, lattice, model.measurements()); // or mps.measure(...);
      // maquis::measurements results(mps, lattice, model.measurements());

      /// Write results
      storage::archive ar(parms["resultfile"].str(), "w");
      ar["/spectrum/results"] << results;
      ar["/parameters"] << parms;
      ar["/parameters"] << model_parms;

    } catch (std::exception & e) {
      maquis::cerr << "Exception thrown!" << std::endl;
      maquis::cerr << e.what() << std::endl;
      return -1;
  }
}



template <class Matrix, class Symm>
class mps
{
public:
  mps(parameters p);
  results_type measure();
};


template <class Matrix>
class mps_base<Matrix>
{
public:
  virtual results_type measure()=0;
};

template <class Matrix, class Symm>
class mps_model : public mps_base<Matrix>
{
public:
  mps_model(mps<Matrix,Symm> const& m) : model(p) {}
  results_type measure() { return model.measure();}
private:
  mps<Matrix,Symm> model;
};

template <class Matrix>
class mps<Matrix,DefaultSymm> 
{
public:
  mps(parameters p)
  {
    switch(....)
    {
      case '0': impl = boost::shared_ptr<mps_model<Matrix,NoSymm> >(mps_model<Matrix,NoSymm>(p));
    }
  }
  results_type measure() { return impl->measure();}
private:
  boost::shared_ptr<mps_base> impl;
};



////////




template <class Matrix, class Symm>
class mps : public mps_base<Matrix>
{
public:
  mps(parameters p);
  results_type measure();
};


template <class Matrix>
class mps_base<Matrix>
{
public:
  virtual results_type measure()=0;
};

template <class Matrix>
class mps<Matrix,DefaultSymm>
{
public:
  mps(parameters p)
  {
    switch(....)
    {
      case '0': impl = boost::shared_ptr<mps<Matrix,NoSymm> >(mps<Matrix,NoSymm>(p));
    }
  }
  results_type measure() { return impl->measure();}
private:
  boost::shared_ptr<mps_base> impl;
};





