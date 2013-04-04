#ifndef AMBIENT_UTILS_MEMORY
#define AMBIENT_UTILS_MEMORY
#include "ambient/utils/singleton.hpp"
#include "boost/pool/singleton_pool.hpp"

namespace ambient { namespace memory {

    template<size_t S>
    class region {
    public:
        region(){
            this->buffers.push_back(std::malloc(S));
            this->buffer = &this->buffers[0];
            this->iterator = (char*)*this->buffer;
        }
        void realloc(){
            if(*this->buffer == this->buffers.back()){
                this->buffers.push_back(std::malloc(S));
                this->buffer = &this->buffers.back();
            }else
                this->buffer++;
            this->iterator = (char*)*this->buffer;
        }
        void* malloc(size_t sz){
            if(((size_t)iterator + sz - (size_t)*this->buffer) >= S) realloc();
            void* m = (void*)iterator;
            iterator += sz;
            return m;
        }
        void reset(){
            this->buffer = &this->buffers[0];
            this->iterator = (char*)*this->buffer;
        }
    private:
        std::vector<void*> buffers;
        void** buffer;
        char* iterator;
    };

    class bulk : public singleton< bulk >
    {
        // for shared heap: (slightly improves consumption)
        // typedef boost::details::pool::default_mutex mutex;
        // boost::details::pool::guard<mutex> g(mtx);
    public:
        bulk(){
            this->arity = ambient::get_num_threads();
            this->set   = new region<AMBIENT_BULK_CHUNK>[arity];
        }
       ~bulk(){
            delete[] this->set;
        }
        template<size_t S>
        void* malloc(){
            return set[AMBIENT_THREAD_ID].malloc(S);
        }
        void* malloc(size_t sz){
            return set[AMBIENT_THREAD_ID].malloc(sz);
        }
        void drop(){
            for(int i = 0; i < arity; i++) this->set[i].reset();
        }
    private:
        region<AMBIENT_BULK_CHUNK>* set;
        int arity;
    };

} }

namespace ambient {
    extern memory::bulk& bulk;
}

namespace ambient { namespace memory {

    class pool : public singleton< pool >
    {
    private:
        struct empty {};
    public:
        pool(){
            this->arity = ambient::get_num_threads();
        }
       ~pool(){
        }
        void* malloc(size_t sz, int r){
            if(r == 0 && sz < AMBIENT_BULK_CHUNK) return ambient::bulk.malloc(sz);
            return std::malloc(sz); 
        }
        void free(void* ptr, size_t sz, int r){
            if(ptr == NULL || (r == 0 && sz < AMBIENT_BULK_CHUNK)) return;
            return std::free(ptr);
        }
        template<size_t S>
        static void* malloc(){
             return boost::singleton_pool<empty, S, boost::default_user_allocator_new_delete, boost::details::pool::null_mutex>::malloc();
        }
        template<size_t S>
        static void free(void* ptr){
             boost::singleton_pool<empty, S, boost::default_user_allocator_new_delete, boost::details::pool::null_mutex>::free(ptr);
        }
        template<typename T>
        static void* malloc(){
             return malloc<sizeof(T)>();
        }
        template<typename T>
        static void free(void* ptr){
             return free<sizeof(T)>(ptr);
        }
    private:
        int arity;
    };

} }

namespace ambient {
    extern memory::pool& pool;
}

#endif
