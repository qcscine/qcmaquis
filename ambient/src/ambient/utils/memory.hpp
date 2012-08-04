#ifndef AMBIENT_UTILS_MEMORY
#define AMBIENT_UTILS_MEMORY
#include "ambient/utils/singleton.hpp"
#include "boost/pool/singleton_pool.hpp"

#define BULK_LENGTH 8388608*16

namespace ambient { namespace utils {

    struct empty { };

    class bulk_memory : public singleton< bulk_memory > 
    {
    public:
        bulk_memory(){
            this->pool = malloc(BULK_LENGTH);
            this->iterator = (char*)this->pool;
        }
       ~bulk_memory(){
            free(this->pool);
        }
        template<size_t S>
        void* get(){
            void* result = this->iterator;
            this->iterator += 16*((size_t)(S/16)+1);
            return result;
        }
        void refresh(){
            this->iterator = (char*)this->pool;
        }
    public:
        char* iterator;
        void* pool;
    };

    template<typename T>
    inline void intrusive_ptr_add_ref(T* p){
        ++(p->references);
    }

    template<typename T>
    inline void intrusive_ptr_release(T* p){
        if(--(p->references) == 0) delete p;
    } 

    template <size_t size>
    class container {
    public:
        inline container()
        : references(0){ }
    private:
        char memory[size];
        long references;
        friend void intrusive_ptr_add_ref<>(container* p);
        friend void intrusive_ptr_release<>(container* p);
    };

} }

namespace ambient {
    extern utils::bulk_memory& bulk_pool;
}

#endif
