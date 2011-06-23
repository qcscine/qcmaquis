
namespace align{

#define MEM_BOUNDARY 16 // 2**7 = 128, align on 128 bits

#if (defined(__ICL) || defined(_MSC_VER) || defined(__ICC))
#include <fvec.h>
    inline void *aligned_malloc (size_t size, size_t align = 16)  {  return _mm_malloc(size+align,align);  }
    inline void  aligned_free   (void *p)                    {  return _mm_free(p); }
#elif defined (__CYGWIN__)
#include <malloc.h>
    inline void *aligned_malloc (size_t size, size_t align = 16)  {  return __mingw_aligned_malloc(size+align,align);  }
    inline void  aligned_free   (void *p)                    {  return __mingw_aligned_free(p);             }
#elif defined(__FreeBSD__)
#include <stdlib.h>
    inline void* aligned_malloc (size_t size, size_t align = 16) {  return malloc(size); }
    inline void  aligned_free   (void *p)                   {  return free(p); } 
#elif defined(__APPLE__)
#include <xmmintrin.h>
    inline void* aligned_malloc (size_t size, size_t align = 16) {  return _mm_malloc(size+align,align); }
    inline void  aligned_free   (void *p)                   {  return _mm_free(p); }
#else 
#include <malloc.h>
    inline void* aligned_malloc (size_t size, size_t align = 16) {  return memalign(align,size+align); }
    inline void  aligned_free   (void *p)                   {  return free(p); }
#endif    
    
}// end namespace
