//#define TOUCHSTACK_ACCESS_CHECK

#include <cilk/cilk.h>

namespace ambient{

    template<typename T>
    inline touchstack<T>::touchstack(){
        ri = wi = content = (T*)malloc(TOUCHSTACK_LENGTH*sizeof(T));
    }

    template<typename T>
    inline touchstack<T>::~touchstack(){
        free(content);
    }

    template<typename T>
    inline T touchstack<T>::pick(){
        return *ri++;
    }

    template<typename T>
    inline T touchstack<T>::back(){
        return *(wi-1);
    }

    template<typename T>
    inline void touchstack<T>::push_back(T e){
        *wi++ = e;
#ifdef TOUCHSTACK_ACCESS_CHECK
        if(this->size() == TOUCHSTACK_LENGTH){
            printf("\n\n\n\n\n\n\n\nERROR: END OF TOUCHSTACK HAS BEEN REACHED (%d)!\n\n\n\n\n\n\n", TOUCHSTACK_LENGTH);
        }
#endif
    }

    template<typename T>
    inline bool touchstack<T>::end_reached(){
        return (ri == wi);
    }

    template<typename T>
    inline void touchstack<T>::repeat(){
        ri = content;
    }

    template<typename T>
    inline void touchstack<T>::reset(){
        ri = wi = content;
    }

    template<typename T>
    inline bool touchstack<T>::empty(){
        return (wi == content);
    }

    template<typename T>
    inline size_t touchstack<T>::size(){
        return ((size_t)wi-(size_t)content)/sizeof(T);
    }

    template<typename T>
    inline void touchstack<T>::sort(){ // insertion sort
        int length = ((size_t)wi - (size_t)content)/sizeof(T);
        for(int i = 1; i < length; i++){
            T value = content[i];
            int j = i - 1;
            bool done = false;
            do {
                if(content[j]->get_weight() < value->get_weight()){
                    content[j + 1] = content[j];
                    if(--j < 0) done = true;
                }else
                    done = true;
            } while(!done);
            content[j + 1] = value;
        }
    }

}
