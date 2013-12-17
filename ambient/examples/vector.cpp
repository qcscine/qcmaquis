#include "ambient/ambient.hpp"

template<typename T> class vector;

namespace detail { using namespace ambient;

    template<typename T>
    void init_value(unbound< vector<T> >& a, T& value){
        size_t size = get_length(a)-1;
        T* a_ = updated(a).elements;
        for(size_t i = 0; i < size; ++i) a_[i] = value;
        updated(a).count = 0;
    }

    template<typename T>
    void add(vector<T>& a, const vector<T>& b){
        size_t size = get_length(a)-1;
        T* a_ = current(a).elements;
        T* b_ = current(b).elements;
        T* result = updated(a).elements;
        for(size_t i = 0; i < size; ++i) result[i] = a_[i] + b_[i];
        updated(a).count = current(a).count+1;
    }
}

ambient_reg(detail::init_value, init_value)
ambient_reg(detail::add, add)

template <typename T>
class vector {
public:
    ambient_version(
        T count;
        T elements[ AMBIENT_VAR_LENGTH ];
    );
    vector(size_t length, T value) : versioned(length+1, sizeof(T)) {
        init_value<T>::spawn(*this, value);
    }
};

int main(){
    using namespace ambient;

    vector<int> a(10, 13);
    vector<int> b(10, 10);
                                 { scope<> select(1);
    add<int>::spawn(a, b);
                                 }
                                 { scope<> select(0); 
    add<int>::spawn(a, b);
                                 }
    for(int i = 0; i < 10; i++)
    cout << "After sync: " << get(a).elements[i] << "; count: " << get(a).count << "\n";
    return 0;
}

