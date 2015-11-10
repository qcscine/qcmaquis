#include "ambient/ambient.hpp"
#include "ambient/container/future.hpp"

int main(){
    ambient::future<double> sum;
    int a = 10;
    int b = 20;

    ambient::bind_cpu([](ambient::future<double>& sum_, int a, int b){
        sum_.set(a+b);
    }, sum, a, b);

    std::cout << "Future value: " << sum << "\n";
    return 0;
}
