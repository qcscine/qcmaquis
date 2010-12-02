/*
 *  Tim_allocator.h
 *  TEST_ALIGNSTL
 *
 *  Created by Tim Ewart on 03.11.10.
 *  Copyright 2010 University of Geneva. All rights reserved.
 *
 */


#ifndef __ALLOCATOR__
#define __ALLOCATOR__

#if (defined(__ICL) || defined(_MSC_VER) || defined(__ICC))
#include <fvec.h>
inline void *aligned_malloc (size_t size, size_t align)  {  return _mm_malloc(size+align,align);  }
inline void  aligned_free   (void *p)                    {  return _mm_free(p); }
#elif defined (__CYGWIN__)
#include <xmmintrin.h>
inline void *aligned_malloc (size_t size, size_t align)  {  return _mm_malloc(size+align,align);  }
inline void  aligned_free   (void *p)                    {  return _mm_free(p); }
#elif defined(__MINGW32__)
#include <malloc.h>
inline void *aligned_malloc (size_t size, size_t align)  {  return __mingw_aligned_malloc(size+align,align);  }
inline void  aligned_free   (void *p)                    {  return __mingw_aligned_free(p);             }
#elif defined(__FreeBSD__)
#include <stdlib.h>
inline void* aligned_malloc (size_t size, size_t align) {  return malloc(size); }
inline void  aligned_free   (void *p)                   {  return free(p); } 
#elif defined(__APPLE__)
#include <stdlib.h>
//previous version
//inline void* aligned_malloc (size_t size, size_t align) {  return _mm_malloc(size+align,align); }
//inline void  aligned_free   (void *p)                   {  return _mm_free(p); }
inline void* aligned_malloc (size_t size, size_t align) {  return malloc(size); }
inline void  aligned_free   (void *p)                   {  return free(p); }
#else 
#include <malloc.h>
inline void* aligned_malloc (size_t size, size_t align) {  return memalign(align,size+align); }
inline void  aligned_free   (void *p)                   {  return free(p); }
#endif


template<class T, std::size_t N=16> class alignment_allocator
{
public:
    typedef T value_type;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
	
    typedef T* pointer;
    typedef const T* const_pointer;
	
    typedef T& reference;
    typedef const T& const_reference;
	
public:
    inline alignment_allocator() throw() {}
    template <class T2> inline alignment_allocator(const alignment_allocator<T2,N>&) throw() {}
	
    inline ~alignment_allocator() throw() {}
	
    inline pointer       address(reference       r)       
	{
		return &r; 
	}
    
	inline const_pointer address(const_reference r) const
	{
		return &r; 
	}
	
    inline pointer allocate(size_type n) 
	{ 
		return reinterpret_cast<pointer>(aligned_malloc(n*sizeof(value_type),N)); 
	}
   
	inline void deallocate(pointer p, size_type)
	{
		aligned_free(p); 
	}
	
    inline void construct (pointer p,const value_type& wert)  
	{ 
		new (p) value_type(wert); 
	}
    
	inline void destroy   (pointer p)  
	{	
		p->~value_type();         
	}
	
    inline size_type max_size() const throw() 
	{ 
		return size_type(-1)/sizeof(value_type); 
	}
	
    template<class T2> struct rebind 
	{
		typedef alignment_allocator<T2,N> other; 
	};
};

#endif