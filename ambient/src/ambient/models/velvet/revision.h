#ifndef AMBIENT_MODELS_VELVET_REVISION
#define AMBIENT_MODELS_VELVET_REVISION

#include <atomic>

namespace ambient { namespace models { namespace velvet {

    class revision
    {
    public:
        revision(size_t extent, void* g, ambient::locality l);
        void* operator new (size_t size);
        void  operator delete (void* ptr);
        template<typename T> operator T* (){ return (T*)data; }

        void embed(void* memory);
        void reuse(revision& r);

        void use();
        void release();
        void complete();

        bool locked();
        bool valid();

        size_t    extent;
        revision* parent;
        void*     generator;
        void*     transfer;
        void*     data;
        int       sid;
        std::atomic<int> users;
        int       region;
        ambient::locality state;
    };

} } }

#endif
