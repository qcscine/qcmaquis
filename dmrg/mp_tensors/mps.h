#ifndef MPS_H
#define MPS_H

#include "mp_tensors/mpstensor.h"
#include "mp_tensors/mpotensor.h"

template<class Matrix, class SymmGroup>
class MPS : public std::vector<MPSTensor<Matrix, SymmGroup> >
{
public:
    typedef std::size_t size_t;
    
    MPS(size_t L, size_t Mmax, Index<SymmGroup> phys);
    
    size_t length() const { return this->size(); }
    Index<SymmGroup> const & site_dim(size_t i) const { return (*this)[i].site_dim(); }
    Index<SymmGroup> const & row_dim(size_t i) const { return (*this)[i].row_dim(); }
    Index<SymmGroup> const & col_dim(size_t i) const { return (*this)[i].col_dim(); }
    
    void canonize(size_t center);
    block_matrix<Matrix, SymmGroup> canonize_left_step(size_t site);
    block_matrix<Matrix, SymmGroup> canonize_right_step(size_t site);
    
    void normalize_left();
    void normalize_right();
    
    Boundary<Matrix, SymmGroup> start_mtx() const;
    
    std::string description() const;
    
    void stupid_grow_pair(std::size_t l, double alpha, double cutoff);
    void grow_l2r_sweep(MPOTensor<Matrix, SymmGroup> const & mpo,
                        Boundary<Matrix, SymmGroup> const & left,
                        Boundary<Matrix, SymmGroup> const & right,
                        std::size_t l, double alpha,
                        double cutoff, std::size_t Mmax);
    
private:
    typename Matrix::value_type canonize_left();
    typename Matrix::value_type canonize_right();
    
};

#include "mp_tensors/mps.hpp"

#endif
