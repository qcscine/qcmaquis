#ifndef AMBIENT_CONTROLLERS_VELVET_MEMORY
#define AMBIENT_CONTROLLERS_VELVET_MEMORY
#include "ambient/utils/singleton.hpp"
#include "boost/pool/singleton_pool.hpp"

#define BULK_LENGTH 8388608*224

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
        template<class T>
        void* get(){
            void* result = this->iterator;
            this->iterator += 16*((size_t)(sizeof(T)/16)+1);
            return result;
        }
        void refresh(){
            this->iterator = (char*)this->pool;
        }
    public:
        char* iterator;
        void* pool;
    };

} }

namespace ambient {
    extern utils::bulk_memory& bulk_pool;
}

#endif

