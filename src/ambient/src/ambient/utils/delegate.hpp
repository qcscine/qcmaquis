#ifndef AMBIENT_UTILS_DELEGATE
#define AMBIENT_UTILS_DELEGATE

namespace ambient {

    class delegate {
    public:
        inline delegate() 
        : length(0), callbacks(NULL) 
        {
        }

        template<typename T>
        inline void operator+=(void(*callback)(T&)){
            callbacks = (void(**)())realloc(this->callbacks, (++this->length)*sizeof(void(*)()));
            callbacks[length-1] = (void(*)())callback;
        }

        template<typename T>
        inline void operator()(T& arg){
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
