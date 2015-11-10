#include "ambient/ambient.hpp"
#include "ambient/container/atomic.hpp"
#include "ambient/utils/reduce.hpp"

int main(){
    const int N = 100;
    std::vector<ambient::atomic<int> > a(N, ambient::atomic<int>(1));

    ambient::reduce(a, [](ambient::atomic<int>& a, const ambient::atomic<int>& b){
        ambient::bind_cpu([](ambient::atomic<int>& dst, const ambient::atomic<int>& src){
            dst.set(dst.get()+src.get());
        }, a, b);
    });

    std::cout << "Reduced value: " << ambient::load(a[0]).get() << "\n";
    return 0;
}
