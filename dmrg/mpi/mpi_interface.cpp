//
// Copyright (C) 2019 Quantum Simulation Technologies, Inc. - All Rights Reserved
//

#include <iostream>
#include <iomanip>
#include <cassert>
#include <thread>
#include <stdexcept>
#include <array>
#include "mpi_interface.h"

using namespace std;
using namespace maquis;

mpiInterface::mpiInterface()
 : depth_(0), cnt_(0), nprow_(0), npcol_(0), context_(0), myprow_(0), mypcol_(0) {

#ifdef HAVE_MPI_H
  int provided;
  MPI_Init_thread(0, 0, MPI_THREAD_MULTIPLE, &provided);
  if (provided != MPI_THREAD_MULTIPLE)
    throw runtime_error("MPI_THREAD_MULTIPLE not provided");

  MPI_Comm_rank(MPI_COMM_WORLD, &wrank_);
  MPI_Comm_size(MPI_COMM_WORLD, &wsize_);
  rank_ = wrank_;
  size_ = wsize_;
  // print out the node name
  {
    constexpr const size_t maxlen = MPI_MAX_PROCESSOR_NAME;
    int len;
    char name[maxlen];
    MPI_Get_processor_name(name, &len);

    unique_ptr<char[]> buf(new char[maxlen*size_]);
    unique_ptr<int[]> lens(new int[size_]);
    MPI_Gather(static_cast<void*>(name), maxlen, MPI_CHAR, static_cast<void*>(buf.get()), maxlen, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Gather(static_cast<void*>(&len),      1, MPI_INT,  static_cast<void*>(lens.get()),     1, MPI_INT,  0, MPI_COMM_WORLD);
    if (rank() == 0) {
      for (int i = 0; i != size_; ++i)
        cout << "    " << string(&buf[i*maxlen], &buf[i*maxlen+lens[i]]) << endl;
      cout << endl;
    }
  }

  // set MPI_COMM_WORLD to mpi_comm_
  mpi_comm_ = MPI_COMM_WORLD;
#else
  wrank_ = 0;
  wsize_ = 1;
  rank_ = wrank_;
  size_ = wsize_;
#endif
}


mpiInterface::~mpiInterface() {
#ifdef HAVE_MPI_H
  MPI_Finalize();
#endif
}


void mpiInterface::barrier() const {
#ifdef HAVE_MPI_H
  MPI_Barrier(mpi_comm_);
#endif
}


void mpiInterface::allreduce(double* a, const size_t size) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_DOUBLE, MPI_SUM, mpi_comm_);
#endif
}


void mpiInterface::allreduce(int* a, const size_t size) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_INT, MPI_SUM, mpi_comm_);
#endif
}


void mpiInterface::allreduce(complex<double>* a, const size_t size) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_);
#endif
}


void mpiInterface::broadcast(size_t* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  static_assert(sizeof(size_t) == sizeof(unsigned long long), "size_t is assumed to be the same size as unsigned long long");
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Bcast(static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_UNSIGNED_LONG_LONG, root, mpi_comm_);
#endif
}


void mpiInterface::broadcast(double* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Bcast(static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_DOUBLE, root, mpi_comm_);
#endif
}


void mpiInterface::broadcast(complex<double>* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Bcast(static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_CXX_DOUBLE_COMPLEX, root, mpi_comm_);
#endif
}


void mpiInterface::allgather(const double* send, const size_t ssize, double* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_DOUBLE, static_cast<void*>(rec), rsize, MPI_DOUBLE, mpi_comm_);
#else
  assert(ssize == rsize);
//copy_n(send, ssize, rec);
#endif
}


void mpiInterface::allgather(const complex<double>* send, const size_t ssize, complex<double>* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_CXX_DOUBLE_COMPLEX, static_cast<void*>(rec), rsize, MPI_CXX_DOUBLE_COMPLEX, mpi_comm_);
#else
  assert(ssize == rsize);
//copy_n(send, ssize, rec);
#endif
}


void mpiInterface::allgather(const size_t* send, const size_t ssize, size_t* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  static_assert(sizeof(size_t) == sizeof(unsigned long long), "size_t is assumed to be the same size as unsigned long long");
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_UNSIGNED_LONG_LONG, static_cast<void*>(rec), rsize, MPI_UNSIGNED_LONG_LONG, mpi_comm_);
#else
  assert(ssize == rsize);
//copy_n(send, ssize, rec);
#endif
}


void mpiInterface::allgather(const int* send, const size_t ssize, int* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_INT, static_cast<void*>(rec), rsize, MPI_INT, mpi_comm_);
#else
  assert(ssize == rsize);
//copy_n(send, ssize, rec);
#endif
}


// MPI Communicators
void mpiInterface::split(const int n) {
#ifdef HAVE_MPI_H
  MPI_Comm new_comm;
  const int icomm = rank_ % n;
  mpi_comm_old_.push_back(pair<MPI_Comm,array<int,5>>(mpi_comm_, {context_, nprow_, npcol_, myprow_, mypcol_}));

  ++depth_;

  MPI_Comm_split(mpi_comm_, icomm, wrank_, &new_comm);
  mpi_comm_ = new_comm;
  MPI_Comm_rank(mpi_comm_, &rank_);
  MPI_Comm_size(mpi_comm_, &size_);
#endif
}


void mpiInterface::merge() {
#ifdef HAVE_MPI_H
  MPI_Comm_free(&mpi_comm_);

  --depth_;

  mpi_comm_ = get<0>(mpi_comm_old_[depth_]);
  MPI_Comm_rank(mpi_comm_, &rank_);
  MPI_Comm_size(mpi_comm_, &size_);

  mpi_comm_old_.pop_back();
#endif
}


int mpiInterface::pnum(const int prow, const int pcol) const {
  return 0;
}
