#ifndef AMBIENT_UTILS_MEMORY
#define AMBIENT_UTILS_MEMORY
#include "ambient/utils/singleton.hpp"
#include "boost/pool/singleton_pool.hpp"

#define BULK_LENGTH 8388608*50

namespace ambient { namespace utils {

    struct empty { };

    class bulk_memory : public singleton< bulk_memory > 
    {
    public:
        bulk_memory(){
            this->arity = __cilkrts_get_nworkers();
            this->pools = (void**)malloc(sizeof(void*)*this->arity);
            for(int i = 0; i < this->arity; i++)
                this->pools[i] = malloc(BULK_LENGTH);
            this->iterators = (char**)malloc(sizeof(char*)*this->arity);
            for(int i = 0; i < this->arity; i++)
                this->iterators[i] = (char*)this->pools[i];
        }
       ~bulk_memory(){
            for(int i = 0; i < this->arity; i++)
                free(this->pools[i]);
            free(this->iterators);
            free(this->pools);
        }
        template<size_t S>
        void* get(){
            char*& iterator = this->iterators[__cilkrts_get_worker_number()];
            void* result = iterator;
            iterator += S; // 16*((size_t)(S/16)+1); // alignment variant
            return result;
        }
        void refresh(){
            for(int i = 0; i < this->arity; i++)
                this->iterators[i] = (char*)this->pools[i];
        }
    public:
        char** iterators;
        void** pools;
        int arity;
    };

} }

namespace ambient {
    extern utils::bulk_memory& bulk_pool;
}

#endif
