#ifndef AMBIENT_UTILS_TOUCHSTACK
#define AMBIENT_UTILS_TOUCHSTACK
#define TOUCHSTACK_LENGTH 256*sizeof(T)

namespace ambient{

    template<typename T>
    class touchstack {
    public:
        inline touchstack();
        inline T pick();
        inline void push_back(T e);
        inline bool end_reached();
        inline void reset();
        inline bool empty();
        inline void sort();
    private:
        T content[TOUCHSTACK_LENGTH];
        T* wi; 
        T* ri;
    };

}

#include "ambient/utils/touchstack.hpp"
#endif
