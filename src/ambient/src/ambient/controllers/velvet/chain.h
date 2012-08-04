#ifndef AMBIENT_CONTROLLERS_VELVET_CHAIN
#define AMBIENT_CONTROLLERS_VELVET_CHAIN

namespace ambient { namespace controllers { namespace velvet {

    using ambient::models::velvet::revision;

    class chain {
    public:
        inline void* operator new (size_t size);
        inline void operator delete (void* ptr); 
        inline chain(cfunctor* f);
        inline void push_back(cfunctor* f);
        inline void execute();
        inline bool ready(); 
        inline bool constrains(cfunctor* op);
        std::vector<cfunctor*> content;
    };

} } }

#endif
