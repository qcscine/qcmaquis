#ifndef MPO_H
#define MPO_H

#include <vector>

#include "mp_tensors/mpotensor.h"

template<class Matrix, class SymmGroup>
class MPO : public std::vector<MPOTensor<Matrix, SymmGroup> >
{
public:
    typedef MPOTensor<Matrix, SymmGroup> elem_type;
    
    MPO(std::size_t L, elem_type elem = elem_type())
    : std::vector<elem_type>(L, elem)
    { }
    
    std::size_t length() const { return this->size(); }
};

#endif
