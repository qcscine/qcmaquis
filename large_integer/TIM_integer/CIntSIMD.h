/*
 *  CIntSIMD.h
 *  512BITS
 *
 *  Created by Tim Ewart on 02.11.10.
 *  Copyright 2010 University of Geneva. All rights reserved.
 *
 */

extern const long long int BASE = 0x3FFFFFFF;

__m128i xmm0,xmm1,xmm2;


template <class Type>
class CContainerSIMD
{
public:
	
	
	CContainerSIMD(Type nInit = 0)
	{
		m_TabnData[0]= static_cast<unsigned long long int> (nInit);
		m_TabnData[1]= static_cast<unsigned long long int> (nInit);
	}
	
	CContainerSIMD(const CContainerSIMD<Type>& A_int128)
	{
		m_TabnData[0]= A_int128.m_TabnData[0];
		m_TabnData[1]= A_int128.m_TabnData[1];
	}
	
	long long int &operator[](size_t nIndex)
	{
		return m_TabnData[nIndex];
	}
	
	
	CContainerSIMD<Type>& operator=(const Type A)
	{
		m_TabnData[0] = static_cast<unsigned long long int> (A);
		return *this;
	}
	
	CContainerSIMD<Type>& operator=(const CContainerSIMD<Type>   & A)
	{
		m_TabnData[0] = A.m_TabnData[0];
		m_TabnData[0] = A.m_TabnData[0];
		return *this;
	}
	
	friend CContainerSIMD<Type> operator+(const CContainerSIMD<Type> & A_int128, const CContainerSIMD<Type> & B_int128)
	{
		
		
		CContainerSIMD<Type> C_int128;
		
		xmm0 = _mm_load_si128((__m128i*)(&A_int128.m_TabnData[0]));		// Load A
		xmm1 = _mm_load_si128((__m128i*)(&B_int128.m_TabnData[0]));		// Load B
		
		xmm0 = _mm_add_epi32(xmm0, xmm1);								// Add A +B
		xmm2 = _mm_load_si128(&xmm0);									// Copy of the sum 6 cycles
		
		xmm2 = _mm_srai_epi32(xmm2, 30);								// Get carry bit (>> C operator)
		xmm2 = _mm_srli_si128(xmm2, 8);									// Move the carry bit 
		xmm0 = _mm_and_si128(xmm0,_mm_set1_epi64((__m64)BASE));	    	// Put the carry bit to 0 (mask, & operator with the base defintion) 
		xmm0 = _mm_add_epi32(xmm0,xmm2);								// Final add
		
		_mm_store_si128((__m128i*)(&C_int128.m_TabnData[0]),xmm0);					// Set the final results
		
		return C_int128;
	}	
	
	friend CContainerSIMD<Type> operator+(const CContainerSIMD<Type> & A_int128, const  Type A)
	{
		CContainerSIMD<Type> C_int128;
		
		xmm0 = _mm_load_si128((__m128i*)(&A_int128.m_TabnData[0]));				// Load A
		xmm1 = _mm_load_si128((__m128i*)(_mm_set1_epi64((__m64)A)));	// Load B
		
		xmm0 = _mm_add_epi32(xmm0, xmm1);								// Add A +B
		xmm2 = _mm_load_si128(&xmm0);									// Copy of the sum
		xmm2 = _mm_srai_epi32(xmm2, 30);								// Get carry bit (>> C operator)
		xmm2 = _mm_srli_si128(xmm2, 8);									// Move the carry bit 
		xmm0 = _mm_and_si128(xmm0,_mm_set1_epi64((__m64)BASE));	    	// Put the carry bit to 0 (mask, & operator with the base defintion) 
		xmm0 = _mm_add_epi32(xmm0,xmm2);								// Final add
		
		_mm_store_si128((__m128i*)(&C_int128.m_TabnData[0]),xmm0);					// Set the final results
		
		return C_int128;
		
	}
	
	friend CContainerSIMD<Type> operator+(const Type A, const CContainerSIMD<Type>& B_int128)
	{
		CContainerSIMD<Type> C_int128;
		
		xmm0 = _mm_load_si128((__m128i*)(&B_int128.m_TabnData[0]));				// Load A
		xmm1 = _mm_load_si128((__m128i*)(_mm_set1_epi64((__m64)A)));	// Load B
		
		xmm0 = _mm_add_epi32(xmm0, xmm1);								// Add A +B
		xmm2 = _mm_load_si128(&xmm0);									// Copy of the sum
		xmm2 = _mm_srai_epi32(xmm2, 30);								// Get carry bit (>> C operator)
		xmm2 = _mm_srli_si128(xmm2, 8);									// Move the carry bit 
		xmm0 = _mm_and_si128(xmm0,_mm_set1_epi64((__m64)BASE));	    	// Put the carry bit to 0 (mask, & operator with the base defintion) 
		xmm0 = _mm_add_epi32(xmm0,xmm2);								// Final add
		
		_mm_store_si128((__m128i*)(&C_int128.m_TabnData[0]),xmm0);					// Set the final results
		
		return C_int128;
	}
	
	private	:
	
	long long int __attribute__((aligned(16))) m_TabnData[2]; // align in the stack, faster ? public to keep const in overloding of operator +
	
};




