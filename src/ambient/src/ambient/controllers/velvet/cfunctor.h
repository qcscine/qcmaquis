#ifndef AMBIENT_CONTROLLERS_VELVET_CFUNCTOR
#define AMBIENT_CONTROLLERS_VELVET_CFUNCTOR

namespace ambient { namespace controllers { namespace velvet {

    class cfunctor : public models::velvet::sfunctor {
    public:
        inline void* operator new (size_t size);
        inline void operator delete (void* ptr);
        virtual void weight()      = 0;
        virtual void logistics()   = 0;
        virtual void computation() = 0;
        virtual ~cfunctor();
        inline cfunctor();
        inline size_t get_weight(){ return 0; } //credit; }
        inline void set_weight(size_t c){ } //credit = c; }
        inline bool ready();
    };

} } }

#endif
