/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef ITERATOR_BLAS1_H
#define ITERATOR_BLAS1_H

#include <boost/numeric/bindings/blas.hpp>
#include <algorithm>

template<class ForwardIterator, class T>
void iterator_axpy(ForwardIterator in1, ForwardIterator in2,
                   ForwardIterator out1, T val)
{
    std::transform(in1, in2, out1, std::bind2nd(std::multiplies<T>(), val));
}

inline void iterator_axpy(double const * in1, double const * in2,
                          double * out1, double val)
{
    fortran_int_t one = 1, diff = in2-in1;
    daxpy_(&diff, &val, in1, &one, out1, &one);
}

inline void iterator_axpy(std::complex<double> const * in1, std::complex<double> const * in2,
                          std::complex<double> * out1, double val)
{
    throw std::runtime_error("Not implemented.");
}

#endif
