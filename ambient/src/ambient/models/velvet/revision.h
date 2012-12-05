#ifndef AMBIENT_MODELS_VELVET_REVISION
#define AMBIENT_MODELS_VELVET_REVISION

//#include <stdatomic.h>

namespace ambient { namespace models { namespace velvet {

    class revision
    {
    public:
        revision(memspec*, void* g);
        void* operator new (size_t size);
        void operator delete (void* ptr);
        operator char* (){ return (char*)this->data; }
        operator double* (){ return (double*)this->data; }
        operator std::complex<double>* (){ return (std::complex<double>*)this->data; }
        void embed(void* memory, size_t bound);
        void swap(revision&);
        void* get_memory();
        bool valid();
        bool occupied();
        void release();
        void use();

        revision* get_parent(){ return parent; }
        void* get_generator();
        void reset_generator();

        int users; // std::atomic<int>
        size_t sid;
        memspec* spec;
        void* header;
        void* data;
        bool clean;
        void* generator;
        revision* parent;
    };

} } }

#endif
