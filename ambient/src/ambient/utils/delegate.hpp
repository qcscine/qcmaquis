#ifndef AMBIENT_UTILS_DELEGATE
#define AMBIENT_UTILS_DELEGATE
#include <stdlib.h>

namespace ambient {

    class delegate {
    public:
        delegate() 
        : length(0), callbacks(NULL) 
        {
        }

        template<typename T>
        void operator+=(void(*callback)(T&)){
            callbacks = (void(**)())realloc(this->callbacks, (++this->length)*sizeof(void(*)()));
            callbacks[length-1] = (void(*)())callback;
        }

        template<typename T>
        void operator()(T& arg){
            for(size_t i = 0; i < this->length; i++){
                ((void(*)(T&))this->callbacks[i])(arg);
            }
        }

        size_t length;
    private:
        void(**callbacks)();
    };

}

#endif
