#include <iostream>
#include <vector>
#include <emmintrin.h>

#include  <ctime>

extern const long long int BASE = 0x3FFFFFFF;

__m128i xmm0,xmm1,xmm2;

template <class Type = std::size_t, std::size_t NBit = 128>
class Clint
{
public:
	
	Clint(Type NInit = 0):nNBit_(NBit)
	{
		TabnData_[0]=static_cast<unsigned long long int> (NInit);
		TabnData_[1]=0;
	}
	
	Clint(const Clint<Type,NBit>& A_int128 ):nNBit_(NBit)
	{
		nNBit_=A_int128.GetSizeBits();
		TabnData_[0]=A_int128.TabnData_[0];
		TabnData_[1]=A_int128.TabnData_[1];
	}

	~Clint()
	{
	}
	
	
	size_t GetSizeBits() const
	{
		return nNBit_;
	}	
	
	long long int &operator[](size_t nIndex)
	{
		return TabnData_[nIndex];
	}
	
	Clint<Type, NBit>& operator=(Type A)
	{
		TabnData_[0] = static_cast<unsigned long long int> (A);
		return *this;
	}

	long long int __attribute__((aligned(16))) TabnData_[2]; // align in the stack, faster ? public to keep const in overloding of operator +
	
private:

	size_t nNBit_; 

};


template <class Type, std::size_t NBit>
Clint<Type, NBit> operator+(const Clint<Type, NBit> & A_int128, const Clint<Type, NBit> & B_int128)
{

	
	Clint<Type, NBit> C_int128;
	
	xmm0 = _mm_load_si128((__m128i*)(&A_int128.TabnData_[0]));		// Load A
	xmm1 = _mm_load_si128((__m128i*)(&B_int128.TabnData_[0]));		// Load B
	
	xmm0 = _mm_add_epi32(xmm0, xmm1);								// Add A +B
	xmm2 = _mm_load_si128(&xmm0);									// Copy of the sum 6 cycles
	
	xmm2 = _mm_srai_epi32(xmm2, 30);								// Get carry bit (>> C operator)
	xmm2 = _mm_srli_si128(xmm2, 8);									// Move the carry bit 
	xmm0 = _mm_and_si128(xmm0,_mm_set1_epi64((__m64)BASE));	    	// Put the carry bit to 0 (mask, & operator with the base defintion) 
	xmm0 = _mm_add_epi32(xmm0,xmm2);								// Final add
	
	_mm_store_si128((__m128i*)(&C_int128.TabnData_[0]),xmm0);					// Set the final results
	
	return C_int128;
}

template <class Type, std::size_t NBit>
Clint<Type, NBit> operator+(Clint<Type,NBit>& A_int128, Type A)
{
	Clint<> C_int128;

	xmm0 = _mm_load_si128((__m128i*)(&A_int128[0]));				// Load A
	xmm1 = _mm_load_si128((__m128i*)(_mm_set1_epi64((__m64)A)));	// Load B
 
	xmm0 = _mm_add_epi32(xmm0, xmm1);								// Add A +B
	xmm2 = _mm_load_si128(&xmm0);									// Copy of the sum
    xmm2 = _mm_srai_epi32(xmm2, 30);								// Get carry bit (>> C operator)
    xmm2 = _mm_srli_si128(xmm2, 8);									// Move the carry bit 
    xmm0 = _mm_and_si128(xmm0,_mm_set1_epi64((__m64)BASE));	    	// Put the carry bit to 0 (mask, & operator with the base defintion) 
	xmm0 = _mm_add_epi32(xmm0,xmm2);								// Final add
 
    _mm_store_si128((__m128i*)(&C_int128[0]),xmm0);					// Set the final results
 
    return C_int128;
 
}


template <class Type, std::size_t NBit>
Clint<Type, NBit> operator+(Type A, Clint<Type,NBit>& B_int128)
{
	Clint<> C_int128;
	
	xmm0 = _mm_load_si128((__m128i*)(&B_int128[0]));				// Load A
	xmm1 = _mm_load_si128((__m128i*)(_mm_set1_epi64((__m64)A)));	// Load B
	
	xmm0 = _mm_add_epi32(xmm0, xmm1);								// Add A +B
	xmm2 = _mm_load_si128(&xmm0);									// Copy of the sum
    xmm2 = _mm_srai_epi32(xmm2, 30);								// Get carry bit (>> C operator)
    xmm2 = _mm_srli_si128(xmm2, 8);									// Move the carry bit 
    xmm0 = _mm_and_si128(xmm0,_mm_set1_epi64((__m64)BASE));	    	// Put the carry bit to 0 (mask, & operator with the base defintion) 
	xmm0 = _mm_add_epi32(xmm0,xmm2);								// Final add
	
    _mm_store_si128((__m128i*)(&C_int128[0]),xmm0);					// Set the final results
	
    return C_int128;
	
}
 	


int main (int argc, char * const argv[])
{
	Clint<> A,B,C;

	srand(3);
		
	A[0] = rand()%BASE;
	A[1] = rand()%BASE;
	
	C = A+B;
	
    return 0;
}
