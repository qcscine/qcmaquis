#ifndef __ALPS_VECTOR_INTERFACE_HPP__
#define __ALPS_VECTOR_INTERFACE_HPP__

#include <vector>

namespace blas
{
    // This macro creates free functions that call member functions with the same
    // name, e.g. swap_columns(A,i,j) -> A.swap_columns(i,j)
#define IMPLEMENT_FORWARDING(RET,NAME,ARGS,VARS) \
template <typename Vector> \
RET NAME ARGS \
{ \
return m.NAME VARS; \
} 
    IMPLEMENT_FORWARDING(void, resize, (Vector& m, typename Vector::size_type i1), (i1) )
    
}

#undef IMPLEMENT_FORWARDING

#endif

