#ifndef AMBIENT_CONTROLLERS_VELVET_CFUNCTOR
#define AMBIENT_CONTROLLERS_VELVET_CFUNCTOR

namespace ambient { namespace controllers { namespace velvet {

    class cfunctor : public models::velvet::sfunctor {
    public:
        virtual void invoke()        = 0;
        virtual bool ready()         = 0;
        virtual const char* symbol() = 0;
        void push_back(cfunctor* d);
        std::vector<cfunctor*> deps;
    };

} } }

#endif
