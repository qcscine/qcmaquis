#include <iostream>
#include <iterator>
#include <emmintrin.h>
using std::endl;
using std::cout;

#include "cow_vector.h"
#include "aligned_allocator.h"


// I could've just used a static_cast<...>, but this saves me repeating the type ;-)
template<class T> const T constify(T foo) { return foo; }

class SIC
{
public:
    int val;
    SIC() { cout << "Scream!" << endl; }
    SIC(const SIC& o) { cout << "Copied!" << endl; }
};

int main()
{
    {
        copy_on_write_vector<double> s(100), f;
        f = s;
        
        // this will not make it unique
        double d0 = constify(s)[0];
        cout << "-" << endl;
        // this, however, will!
        double d1 = s[0];
        
        cout << "###" << endl;
        
        // same here...
        cout << constify(s)[0] << endl;
        cout << s[0] << endl;
    }
	
    {
        copy_on_write_vector<double, aligned_allocator<double> > foo(10);
        foo[0] = 1;
        copy_on_write_vector<double, aligned_allocator<double, 8> > bar;
        copy_on_write_vector<double, aligned_allocator<double, 16> > bar2 = bar;
    }
    
    {
        copy_on_write_vector<double, aligned_allocator<double, 8> > bar(10);
        std::copy(bar.begin(), bar.end(), std::ostream_iterator<double>(cout, " ")); cout << endl;
        copy_on_write_vector<double, aligned_allocator<double, 8, true> > bar2(10);
        std::copy(bar2.begin(), bar2.end(), std::ostream_iterator<double>(cout, " ")); cout << endl;
    }
    
    {
        copy_on_write_vector<SIC, aligned_allocator<SIC, 8, false> > bar(10);
        copy_on_write_vector<SIC, aligned_allocator<SIC, 8, true> > bar2(10);
    }
	
    {
		__m128d xmm0,xmm1;
		
		copy_on_write_vector<double, aligned_allocator<double, 16, true> > A(2);
		A[0] = 0.1;
		A[1] = 0.2;
		
		copy_on_write_vector<double, aligned_allocator<double, 16, true> > B(2);		
		B[0] = 0.9;
		B[1] = 0.8;
		//load
		xmm0 = _mm_load_pd(&A[0]);
		xmm1 = _mm_load_pd(&B[0]);
		//add
		xmm0 = _mm_add_pd(xmm0,xmm1);
		//store
		_mm_store_pd(&A[0], xmm0);
		
		std::copy(A.begin(), A.end(), std::ostream_iterator<double>(cout, " ")); cout << endl;
		
    }
	
	
	
	
	
	
	
}

