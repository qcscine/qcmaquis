#ifndef AMBIENT_UTILS_TOUCHSTACK
#define AMBIENT_UTILS_TOUCHSTACK
#define STACK_CONTENT_RESERVATION 10

namespace ambient{

    template<typename T>
    class touchstack {
    public:
        inline touchstack();
        inline ~touchstack();
        inline void push_back(T element);
        inline bool end_reached();
        inline bool alt_end_reached();
        inline void sort();
        inline T pick();
        inline T alt_pick();
        inline T back();
        inline void reset();
        inline void alt_reset();
        inline void clean();
        inline bool empty();
    private:
        T* content;
        size_t write_iterator; 
        size_t read_iterator;
        size_t alt_read_iterator;
        size_t length;
        size_t reserved;
    };

}

#include "ambient/utils/touchstack.hpp"
#endif
