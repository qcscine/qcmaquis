/*********************************************************************************
*
* MPI interface for Scine
* inspired and adapted for Scine from the MPI interface for BAGEL by Toru Shiozaki
*
* Copyright (C) 2019         Stefan Knecht  <stknecht@ethz.ch>
*
**********************************************************************************/
#ifndef MPI_INTERFACE_H
#define MPI_INTERFACE_H

#include <stddef.h>
#include <memory>
#include <complex>
#include <mutex>
#include <vector>
#include <map>

#ifdef HAVE_MPI
#include <mpi.h>
#else
using MPI_Comm = int;
#endif

#ifdef PROG_QCM
namespace maquis
{
#endif

class MPI_interface {
  // protected members are accessible from other members of the same class (or from their "friends"), but also from members of their derived classes
  protected:
    bool initialized_; // is the MPI interface initialized?
    int id_gl_; // global rank
    int id_sm_; // local rank (shared-memory communicator)
    int size_gl_; // size of global communicator aka MPI_COMM_WORLD
    int size_sm_; // size of local (shared-memory) communicator

#ifdef HAVE_MPI
    std::vector<MPI_Comm> mpi_comm_;
    // ctag defines which communicator is returned: 0 = global; 1 = shared-memory; 2 = internode
    MPI_Comm& mpi_comm(const int ctag) { return mpi_comm_.at(ctag); };
#endif

  public:
     // default constructor
     MPI_interface(MPI_Comm* commptr, const int id_GL);
     // destructor
    ~MPI_interface();

    // some simple member functions - definition of function body here
    int id_gl() const { return id_gl_; }
    int size_gl() const { return size_gl_; }
    int id_sm() const { return id_sm_; }
    int size_sm() const { return size_sm_; }

};

// declaration that somewhere in the code there is an MPI_interface object called mpi__.
extern MPI_interface* mpi__;

#ifdef PROG_QCM
}
#endif

#endif

