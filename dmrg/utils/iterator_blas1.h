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
#include <boost/lambda/lambda.hpp>

#include <algorithm>

template<class InputIterator, class OutputIterator, class T>
void iterator_axpy(InputIterator in1, InputIterator in2,
                   OutputIterator out1, T val)
{
    using namespace boost::lambda;
    std::transform(in1, in2, out1, out1, _1*val+_2);
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
