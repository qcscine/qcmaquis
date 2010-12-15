/*
 *  general_matrix_gpu_.h
 *  mps
 *
 *  Created by Tim Ewart on 13.12.10.
 *  Copyright 2010 University of Geneva. All rights reserved.
 *
 */

#include "dense_matrix/gpu/gpu_type_macros.h"

namespace blas
{
	
#define MATRIX_MATRIX_MULTIPLY(T) \
template <typename MemoryBlock> \
const general_matrix<T,MemoryBlock> matrix_matrix_multiply(general_matrix<T,MemoryBlock> const& lhs, general_matrix<T,MemoryBlock> const& rhs) \
{ \
assert( lhs.num_columns() == rhs.num_rows() ); \
gpu::matrix_gpu<T> lhs_gpu(lhs); \
gpu::matrix_gpu<T> rhs_gpu(rhs); \
gpu::matrix_gpu<T> result_gpu(num_rows(lhs), num_columns(rhs)); \
 result_gpu = matrix_matrix_multiply( lhs_gpu, rhs_gpu) ; \
return result_gpu; \
}
IMPLEMENT_FOR_ALL_GPU_TYPES(MATRIX_MATRIX_MULTIPLY)
#undef MATRIX_MATRIX_MULTIPLY

#define MATRIX_SVD(T) \
template <typename MemoryBlock> \
void svd(general_matrix<T, MemoryBlock> & M, general_matrix<T, MemoryBlock> & U, general_matrix<T, MemoryBlock>& V, typename associated_diagonal_matrix<general_matrix<T, MemoryBlock> >::type & S) \
{ \
BOOST_CONCEPT_ASSERT((blas::Matrix<general_matrix<T, MemoryBlock> >)); \
typename general_matrix<T, MemoryBlock>::size_type k = std::min(num_rows(M), num_columns(M)); \
resize(U, num_rows(M), k); \
resize(V, k, num_columns(M)); \
std::vector<T> S_(k); \
std::size_t size1 = M.size1(); \
std::size_t size2 = M.size2(); \
gpu::matrix_gpu<T> U_gpu(U); \
gpu::matrix_gpu<T> V_gpu(V); \
gpu::matrix_gpu<T> M_gpu(M); \
gpu::matrix_gpu<T> S_gpu(size1,size2,0); \
 gpu::svd(M_gpu,U_gpu,V_gpu,S_gpu); \
blas::general_matrix<T>  S_matrix(S_gpu); \
for (std::size_t i=0; i < k; i++) \
	S_[i] = S_matrix(i,0);  \ // wierd due to cula ...
U = U_gpu;\
V = V_gpu;\
} 
IMPLEMENT_FOR_ALL_GPU_TYPES(MATRIX_SVD)
#undef MATRIX_SVD	
	
	
	
	
	
	
	
}
