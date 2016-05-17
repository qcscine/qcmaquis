#include "ambient/ambient.hpp"
#include "ambient/container/vector.hpp"

template<typename T>
void reverse(ambient::vector<T>& vec){
    int start = 0;
    int end = vec.size();
    while((start != end) && (start != --end)){
        std::swap(vec[start],vec[end]); start++;
    }
}

int main(){
    ambient::vector<int> a(100);         // zero initialised vector
    ambient::async(reverse<int>, a);  // reverse vector asynchronously
    ambient::async([](ambient::vector<int>& vec){ reverse(vec); }, a);
    ambient::sync();                     // wait for operations to finish
    return 0;
}
