#ifndef AMBIENT_MODELS_VELVET_REVISION
#define AMBIENT_MODELS_VELVET_REVISION

//#include <stdatomic.h>

namespace ambient { namespace models { namespace velvet {

    class revision
    {
    public:
        revision(memspec*, void* g);
        void* operator new (size_t size);
        void  operator delete (void* ptr);
        template<typename T> operator T* (){ return (T*)data; }

        void embed(void* memory, size_t bound);
        void reuse(revision& r);

        void use();
        void release();
        void complete();

        bool locked();
        bool remote();
        bool origin();
        bool valid();

        memspec*  spec;
        revision* parent;
        void*     generator;
        void*     header;
        void*     data;
        size_t    sid;
        int       users;

        enum { VOID,
               PURE, 
               COPY, 
               WAIT, 
               ORIG } state;
    };

} } }

#endif
