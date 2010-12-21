/*
 *  general_matrix_gpu_.h
 *  mps
 *
 *  Created by Tim Ewart on 13.12.10.
 *  Copyright 2010 University of Geneva. All rights reserved.
 *
 */

#include "gpu_type_macros.h"

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

#define MULTIPLIES_ASSIGN(T) \
template <typename MemoryBlock> \
void multiplies_assign(general_matrix<T,MemoryBlock>& m, T const& t) \
{ \
gpu::matrix_gpu<T> m_gpu(m); \
if( !(m.is_shrinkable()) ) \
{ \
	multiplies_assign(m_gpu, t); \
} \
else \
{ \
	exit(0); \
} \
m = m_gpu; \
}
IMPLEMENT_FOR_ALL_GPU_TYPES(MULTIPLIES_ASSIGN)
#undef MULTIPLIES_ASSIGN	
	
#define PLUS_ASSIGN(T) \
template <typename MemoryBlock> \
void plus_assign(general_matrix<T,MemoryBlock>& this_m, general_matrix<T,MemoryBlock> const& rhs) \
{ \
gpu::matrix_gpu<T> this_gpu(this_m); \
gpu::matrix_gpu<T> rhs_gpu(rhs); \
plus_assign(this_gpu, rhs_gpu); \
this_m = this_gpu; \
} 
IMPLEMENT_FOR_ALL_GPU_TYPES(PLUS_ASSIGN)
#undef PLUS_ASSIGN

#define MINUS_ASSIGN(T) \
template <typename MemoryBlock> \
void minus_assign(general_matrix<T,MemoryBlock>& this_m, general_matrix<T,MemoryBlock> const& rhs) \
{ \
gpu::matrix_gpu<T> this_gpu(this_m); \
gpu::matrix_gpu<T> rhs_gpu(rhs); \
minus_assign(this_gpu, rhs_gpu); \
this_m = this_gpu; \
} 
IMPLEMENT_FOR_ALL_GPU_TYPES(MINUS_ASSIGN)
#undef MINUS_ASSIGN

#define SVD(T) \
template <typename MemoryBlock> \
void svd(general_matrix<T, MemoryBlock>  M, general_matrix<T, MemoryBlock> & U, general_matrix<T, MemoryBlock>& V, typename associated_diagonal_matrix<general_matrix<T, MemoryBlock> >::type & S) \
{ \
BOOST_CONCEPT_ASSERT((blas::Matrix<general_matrix<T, MemoryBlock> >)); \
typename general_matrix<T, MemoryBlock>::size_type k = std::min(num_rows(M), num_columns(M)); \
resize(U, num_rows(M), k); \
resize(V, k, num_columns(M)); \
std::size_t size1 = M.size1(); \
std::size_t size2 = M.size2(); \
std::size_t min = k; \
std::size_t max = std::max(num_rows(M), num_columns(M)); \
gpu::matrix_gpu<T> U_gpu(max,max); \
gpu::matrix_gpu<T> V_gpu(min,min); \
gpu::matrix_gpu<T> M_gpu(M); \
gpu::vector_gpu<T> S_gpu(k); \
gpu::svd(M_gpu,U_gpu,V_gpu,S_gpu); \
U = U_gpu; \
V = V_gpu; \
M = M_gpu; \
blas::vector<T>  Sres(k); \
Sres = S_gpu; \
S = typename associated_diagonal_matrix<general_matrix<T, MemoryBlock> >::type(Sres); \
} 
IMPLEMENT_FOR_ALL_GPU_TYPES(SVD)
#undef SVD	
		
}