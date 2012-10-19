#ifndef AMBIENT_MODELS_VELVET_REVISION
#define AMBIENT_MODELS_VELVET_REVISION

//#include <stdatomic.h>

namespace ambient { namespace models { namespace velvet {

    class revision
    {
    public:
        inline void* operator new (size_t size);
        inline void operator delete (void* ptr);
        inline operator char* (){ return (char*)this->data; }
        inline operator double* (){ return (double*)this->data; }
        inline operator std::complex<double>* (){ return (std::complex<double>*)this->data; }
        inline void embed(void* memory, size_t bound);
        inline void swap(revision&);
        inline void* get_memory();
        inline bool valid();
        inline bool occupied();
        inline void release();
        inline void use();

        inline revision* get_parent(){ return parent; }
        inline revision(memspec*, void* g);
        inline void* get_generator();
        inline void reset_generator();

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
