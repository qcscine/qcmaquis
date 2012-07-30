#ifndef AMBIENT_CONTROLLERS_VELVET_MEMORY
#define AMBIENT_CONTROLLERS_VELVET_MEMORY
#include "ambient/utils/singleton.hpp"

#define INSTRUCTION_LENGTH 8388608
#define CHAIN_CHUNK 48
#define CFUNCTOR_CHUNK 224

#include <queue>

namespace ambient { namespace utils { 

    class chain_memory : public singleton< chain_memory > 
    {
    public:
        chain_memory(){
            this->pool = malloc(INSTRUCTION_LENGTH*CHAIN_CHUNK);
            this->iterator = (char*)this->pool;
        }
       ~chain_memory(){
            free(this->pool);
        }
        void* get(){
            this->iterator += CHAIN_CHUNK;
            return this->iterator;
        }
        void reset(){
            this->iterator = (char*)this->pool;
        }
    public:
        char* iterator;
        void* pool;
    };

    class cfunctor_memory : public singleton< cfunctor_memory > 
    {
    public:
        cfunctor_memory(){
            this->pool = malloc(INSTRUCTION_LENGTH*CFUNCTOR_CHUNK);
            this->iterator = (char*)this->pool;
        }
       ~cfunctor_memory(){
            free(this->pool);
        }
        void* get(){
            this->iterator += CFUNCTOR_CHUNK;
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
    extern utils::chain_memory& chain_pool;
    extern utils::cfunctor_memory& cfunctor_pool;
}

#endif

