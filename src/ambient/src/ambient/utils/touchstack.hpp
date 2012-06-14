//#define TOUCHSTACK_ACCESS_CHECK

namespace ambient{

    template<typename T>
    inline touchstack<T>::touchstack(){
        ri = wi = content = (T*)malloc(TOUCHSTACK_LENGTH);
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
    inline void touchstack<T>::push_back(T e){
        *wi++ = e;
#ifdef TOUCHSTACK_ACCESS_CHECK
        if(((size_t)wi-(size_t)content) == TOUCHSTACK_LENGTH){
            printf("ERROR: end of touchstack has been reached!\n");
        }
#endif
    }

    template<typename T>
    inline bool touchstack<T>::end_reached(){
        return (ri == wi);
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
