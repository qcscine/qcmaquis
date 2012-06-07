#ifndef AMBIENT_UTILS_MEMORY
#define AMBIENT_UTILS_MEMORY

namespace ambient {
//    class memory {
//    public:
//       void* alloc(size_t nbytes);
//       void dealloc(void* p);
//    private:
//    } pool;
//
//    void* memory::alloc(size_t nbytes){
//        return malloc(nbytes);
//    }
//
//    void memory::dealloc(void* p){
//        free(p);
//    }
//}

    template<typename T>
    inline void intrusive_ptr_add_ref(T* p){
        ++(p->references);
    }

    template<typename T>
    inline void intrusive_ptr_release(T* p){
        if(--(p->references) == 0) delete p;
    } 

//inline void* operator new(size_t nbytes, ambient::memory& pool){
//    return pool.alloc(nbytes);
//}

}

#endif
