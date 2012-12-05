#ifndef AMBIENT_UTILS_TOUCHSTACK
#define AMBIENT_UTILS_TOUCHSTACK
#define TOUCHSTACK_LENGTH 16388608

// 8192 -- ff_short
// 65536 -- ff_large
// 131072 -- fermi ladder
// 8388608 -- wide fermi ladder

namespace ambient{

    template<typename T>
    class touchstack {
    public:
        touchstack();
       ~touchstack();
        T pick();
        T back();
        void push_back(T e);
        bool end_reached();
        void reset();
        void repeat();
        bool empty();
        size_t size();
        void sort();
    public:
        T* content;
        T* wi; 
        T* ri;
    };

}

#include "ambient/utils/touchstack.hpp"
#endif
