/*****************************************************************************
*
* MPI interface for QCMaquis
*
* Copyright (C) 2019         Stefan Knecht  <stknecht@ethz.ch>
*
*****************************************************************************/
#ifndef MAQUIS_MPI_INTERFACE_H
#define MAQUIS_MPI_INTERFACE_H

#include <stddef.h>
#include <memory>
#include <complex>
#include <mutex>
#include <vector>
#include <map>
#ifdef HAVE_MPI_H
 #include <mpi.h>
#endif

namespace maquis
{


class mpiInterface {
  protected:
    int wrank_;
    int wsize_;
    int rank_;
    int size_;
    int depth_;
    int cnt_;

    int nprow_;
    int npcol_;
    int context_;
    int myprow_;
    int mypcol_;

#ifdef HAVE_MPI_H
    MPI_Comm mpi_comm_;
    std::vector<std::pair<MPI_Comm,std::array<int,5>>> mpi_comm_old_;
#endif

    // maximum size of the MPI buffer - same as, for example, BAGEL uses
    static constexpr size_t bsize = 100000000LU;

  public:
    mpiInterface();
    ~mpiInterface();

    int wrank() const { return wrank_; }
    int wsize() const { return wsize_; }
    int rank() const { return rank_; }
    int size() const { return size_; }
    int depth() const { return depth_; }
    bool last() const { return rank() == size()-1; }

    // collective functions - MPI-1 standard
    // barrier
    void barrier() const;
    // sum reduce and broadcast to each process
    void allreduce(int*, const size_t size) const;
    void allreduce(double*, const size_t size) const;
    void allreduce(std::complex<double>*, const size_t size) const;
    // broadcast
    void broadcast(size_t*, const size_t size, const int root) const;
    void broadcast(double*, const size_t size, const int root) const;
    void broadcast(std::complex<double>*, const size_t size, const int root) const;
    void allgather(const double* send, const size_t ssize, double* rec, const size_t rsize) const;
    void allgather(const std::complex<double>* send, const size_t ssize, std::complex<double>* rec, const size_t rsize) const;
    void allgather(const size_t* send, const size_t ssize, size_t* rec, const size_t rsize) const;
    void allgather(const int* send, const size_t ssize, int* rec, const size_t rsize) const;

#ifdef HAVE_MPI_H
    // communicators. n is the number of processes per communicator.
    const MPI_Comm& mpi_comm() const { return mpi_comm_; }
#endif
    void split(const int n);
    void merge();


    // scalapack
    int nprow() const { return nprow_; }
    int npcol() const { return npcol_; }
    int context() const { return context_; }
    int myprow() const { return myprow_; }
    int mypcol() const { return mypcol_; }

#ifdef HAVE_MPI_H
    // communicators. n is the number of processes per communicator.
    const MPI_Comm& mpi_comm() const { return mpi_comm_; }
#endif

    int pnum(const int prow, const int pcol) const;
    std::pair<int,int> numroc(const int, const int) const;
    std::vector<int> descinit(const int, const int) const;

};

extern mpiInterface* qcm_mpi__;

}


#endif

