#ifndef AMBIENT_CONTROLLERS_VELVET_CFUNCTOR
#define AMBIENT_CONTROLLERS_VELVET_CFUNCTOR

namespace ambient { namespace controllers { namespace velvet {

    class cfunctor : public models::velvet::sfunctor {
    public:
        virtual void invoke()      = 0;
        virtual bool ready(void*)  = 0;
        virtual bool match(void*)  = 0;
        virtual void tag(void*)    = 0;
    };

} } }

#endif
