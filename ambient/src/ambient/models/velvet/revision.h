#ifndef AMBIENT_MODELS_VELVET_REVISION
#define AMBIENT_MODELS_VELVET_REVISION

namespace ambient { namespace models { namespace velvet {

    class sfunctor;

} } }

namespace ambient { namespace controllers { namespace velvet {

    class cfunctor;

} } }

namespace ambient { namespace models { namespace velvet {

    using ambient::models::velvet::sfunctor;
    using ambient::controllers::velvet::cfunctor;

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

        inline revision* get_parent(){ return parent; }
        inline revision(memspec*, bool clean = false);
        inline void* get_generator();
        inline void set_generator(void*);
        inline void reset_generator();

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
