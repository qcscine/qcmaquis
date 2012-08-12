#ifndef AMBIENT_UTILS_TOUCHSTACK
#define AMBIENT_UTILS_TOUCHSTACK
#define TOUCHSTACK_LENGTH 8388608

// 8192 -- ff_short
// 65536 -- ff_large
// 131072 -- fermi ladder
// 8388608 -- wide fermi ladder

namespace ambient{

    template<typename T>
    class touchstack {
    public:
        inline touchstack();
        inline ~touchstack();
        inline T pick();
        inline void push_back(T e);
        inline bool end_reached();
        inline void reset();
        inline void repeat();
        inline void clear();
        inline void purge();
        inline bool empty();
        inline size_t size();
        inline void sort();
    private:
        T* content;
        T* wi; 
        T* ri;
    };

}

#include "ambient/utils/touchstack.hpp"
#endif
