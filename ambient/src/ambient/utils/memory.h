#ifndef AMBIENT_CONTROLLERS_VELVET_MEMORY
#define AMBIENT_CONTROLLERS_VELVET_MEMORY
#include "ambient/utils/singleton.hpp"

#define CHUNK 72
#define LENGTH 8388608

#define INSTRUCTION_LENGTH 8388608
#define INSTRUCTION_CHUNK 72

#include <queue>

namespace ambient { namespace utils { 

    class instruction_memory : public singleton< instruction_memory > 
    {
    public:
        instruction_memory(){
            this->pool = malloc(INSTRUCTION_LENGTH*INSTRUCTION_CHUNK);
            this->iterator = (char*)this->pool;
        }
       ~instruction_memory(){
            free(this->pool);
        }
        void* get(size_t sz){
            //return malloc(sz);
            this->iterator += sz;
            return this->iterator;
        }
        void reset(){
            this->iterator = (char*)this->pool;
        }
    public:
        char* iterator;
        void* pool;
    };

    class memory : public singleton< memory > 
    {
    public:
        inline memory(){
            //this->pool = malloc(LENGTH*CHUNK);
            //this->iterator = (char*)this->pool;
        }
        inline ~memory(){
            //free(this->pool);
        }
        void* get(size_t sz){
            return malloc(sz);
        }
        void reset(void* ptr){
            free(ptr);
        }
    public:
        char* iterator;
        void* pool;
        std::queue<void*> vacants;
    };

} }

namespace ambient {
    extern utils::memory& pool;
    extern utils::instruction_memory& instruction_pool;
}

#endif

