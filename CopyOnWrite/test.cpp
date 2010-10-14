#include <iostream>
#include <iterator>
using std::endl;
using std::cout;

#include "cow_vector.h"
#include "aligned_allocator.h"

// I could've just used a static_cast<...>, but this saves me repeating the type ;-)
template<class T> const T constify(T foo) { return foo; }

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
    }
    
    {
        copy_on_write_vector<double, aligned_allocator<double, 8> > bar(10);
        std::copy(bar.begin(), bar.end(), std::ostream_iterator<double>(cout, " ")); cout << endl;
        copy_on_write_vector<double, aligned_allocator<double, 8, true> > bar2(10);
        std::copy(bar2.begin(), bar2.end(), std::ostream_iterator<double>(cout, " ")); cout << endl;
    }
}
