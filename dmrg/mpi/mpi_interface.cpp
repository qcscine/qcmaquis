/*********************************************************************************
*
* MPI interface for Scine
* inspired and adapted for Scine from the MPI interface for BAGEL by Toru Shiozaki
*
* Copyright (C) 2019         Stefan Knecht  <stknecht@ethz.ch>
*
*********************************************************************************/

#include <iostream>
#include <iomanip>
#include <cassert>
#include <thread>
#include <stdexcept>
#include <array>
#include <mpi_interface.h>

using namespace std;

// QCMaquis specific
#ifdef PROG_QCM
using namespace maquis;
#endif

#ifndef HAVE_MPI
using MPI_Comm = int;
#endif

// specify default constructor
MPI_interface::MPI_interface(MPI_Comm* commptr, const int id_GL)
 : initialized_(true) {

#ifdef HAVE_MPI
  if (!commptr) {
    MPI_Init(0, 0);

    // setup global communication group
    MPI_Comm_rank(MPI_COMM_WORLD, &id_gl_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_gl_);

    mpi_comm_.push_back(MPI_COMM_WORLD);

    // setup shared-memory (NUMA) communication group
    MPI_Comm comm_local;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_local);
    MPI_Comm_size(comm_local, &size_sm_);
    MPI_Comm_rank(comm_local, &id_sm_);

    mpi_comm_.push_back(comm_local);

    // setup communication group for node masters (id_sm==0)
    //  (and all others, but they will not be used ...)
    MPI_Comm_split(MPI_COMM_WORLD, id_sm_ , 0, &comm_local);

    mpi_comm_.push_back(comm_local);

  } else {
    mpi_comm_.push_back(*commptr);
    id_gl_ = id_GL;
  }

#else
  id_gl_   = 0;
  size_gl_ = 1;
  id_sm_   = id_gl_;
  size_sm_ = size_gl_;
#endif
}

// destructor
MPI_interface::~MPI_interface() {
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}

// member functions
//#ifdef HAVE_MPI
//MPI_Comm& MPI_interface::mpi_comm(const int ctag) {
//return mpi_comm_.at(ctag);
//}
//#endif
