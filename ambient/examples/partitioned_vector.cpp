#include "ambient/ambient.hpp"
#include "ambient/container/atomic.hpp"
#include "ambient/container/vector.hpp"
#include "ambient/container/partitioned_vector.hpp"

int main(){
    ambient::partitioned_vector<ambient::vector<int>, 3> a(10, 13);
    ambient::partitioned_vector<ambient::vector<int>, 4> b(10, 39);
    ambient::partitioned_vector<ambient::vector<int>, 5> c(10, 43);

    ambient::for_each(a.begin(), a.end(), [](int& e){ e++; });
    ambient::for_each(a.begin()+1, a.end()-1, [](int& e){ e++; });
    ambient::fill(a.begin(), a.end(), 169);
    ambient::generate(a.begin(), a.end(), []{ return 159; }); 
    ambient::sequence(a.begin(), a.end());
    ambient::transform(a.begin()+2, a.end(), b.begin()+2, c.begin()+2, std::plus<int>());
    ambient::copy(b.begin(), b.begin()+2, c.begin());
    ambient::transform(c.begin(), c.end(), b.begin(), [](int i){ return ++i; });
    ambient::replace(b.begin(), b.end(), 40, 87);

    ambient::partitioned_vector<ambient::vector<int>, 1024*1024> seq(8*1024*1024);
    ambient::sequence(seq.begin(), seq.end());
    auto res = ambient::reduce(seq.begin(), seq.end(), (double)0.);

    ambient::sort(seq.begin(), seq.end());
    ambient::sort(seq.begin(), seq.end(), [](int a, int b){ return a > b; });

    auto pos = ambient::find(seq.begin(), seq.end(), 8*1024*1024-1);

    ambient::sync();
    return 0;
}
