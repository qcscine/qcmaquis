#include <iostream>
using std::endl;
using std::cout;

#include "cow_vector.h"

template<class T> const T constify(T foo) { return foo; }

int main()
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
