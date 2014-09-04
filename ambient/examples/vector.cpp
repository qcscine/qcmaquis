#include "ambient/ambient.hpp"
#include "ambient/container/vector.hpp"

template<typename T>
void reverse(ambient::vector<T>& vec){
    int start = 0;
    int end = get_length(vec);
    while((start != end) && (start != --end)){
        std::swap(vec[start],vec[end]); start++;
    }
}
AMBIENT_EXPORT_TEMPLATE(reverse, reverse_exported);

int main(){
    ambient::vector<int> a(100);      // zero initialised vector
    reverse_exported<int>(a);         // reverse vector asynchronously
    ambient::async(reverse<int>, a);  // ... or without export
    ambient::async([](ambient::vector<int>& vec){ reverse(vec); }, a);
    ambient::sync();                  // wait for operations to finish
    return 0;
}
