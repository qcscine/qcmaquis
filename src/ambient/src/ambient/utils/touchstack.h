#ifndef AMBIENT_UTILS_TOUCHSTACK
#define AMBIENT_UTILS_TOUCHSTACK
#define STACK_CONTENT_RESERVATION 10

namespace ambient{

    template<typename T>
    class touchstack {
    public:
        touchstack();
       ~touchstack();
        void push_back(T element);
        bool end_reached();
        bool alt_end_reached();
        void sort();
        T pick();
        T alt_pick();
        T back();
        void reset();
        void alt_reset();
        void clean();
        bool empty();
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
