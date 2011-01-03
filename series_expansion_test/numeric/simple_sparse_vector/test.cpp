#include "simple_sparse_vector.hpp"
#include <iostream>

using series_expansion::simple_sparse_vector;

int main()
{
    simple_sparse_vector<float> b;
    b(3) += 4;
    b(5) -= 2.455;

    b.print(std::cout);
    for(simple_sparse_vector<float>::iterator it(b.begin()); it != b.end(); ++it)
        std::cout<<it.get_index()<<":"<<*it<<std::endl;
    return 0;
}
