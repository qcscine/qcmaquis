#ifndef AMBIENT_CONTROLLERS_VELVET_CFUNCTOR
#define AMBIENT_CONTROLLERS_VELVET_CFUNCTOR

namespace ambient { namespace controllers { namespace velvet {
    
    using ambient::models::velvet::revision;
    using ambient::models::velvet::transformable;
    using ambient::channels::mpi::request;

    class cfunctor {
        typedef ambient::memory::bulk::allocator<cfunctor*> allocator;
    public:
        virtual void invoke() = 0;
        virtual bool ready() = 0;
        #ifdef AMBIENT_COMPUTATIONAL_DATAFLOW
        virtual const char* name() = 0;
        size_t id(){ return number; }
        size_t number;
        #endif
        void queue(cfunctor* d);
        std::vector<cfunctor*, allocator> deps;
        void** arguments;
    };

    template<class T, int N = 1> class get {};
    template<class T, int N = 0> class set {};

    // {{{ revision get/set
    template<int N>
    class get<revision, N> : public cfunctor {
    public:
        void* operator new (size_t size){ return ambient::bulk.malloc<sizeof(get)>(); }
        void operator delete (void* ptr){ }
        static void spawn(revision& r);
        get(revision& r);
        virtual bool ready();
        virtual void invoke();
        #ifdef AMBIENT_COMPUTATIONAL_DATAFLOW
        virtual const char* name(){ return "get"; }
        #endif
    private:
        revision* target;
        request* handle;
    };
    template<int N>
    class set<revision, N> : public cfunctor {
    public:
        template<class T, int NE> friend class set;
        void* operator new (size_t size, void* placement){ return placement; }
        void* operator new (size_t size){ return ambient::bulk.malloc<sizeof(set)>(); }
        void operator delete (void*, void*){ }
        void operator delete (void* ptr){ }

        static set<revision>& spawn(revision& r);

        set(revision& r);
        template<int NE> set(set<revision,NE>* s);
        virtual void operator >> (int p);
        virtual bool ready();
        virtual void invoke();
        #ifdef AMBIENT_COMPUTATIONAL_DATAFLOW
        virtual const char* name(){ return "set"; }
        #endif
    private:
        bool evaluated;
        std::vector<bool>* states;
        std::vector<request*>* handles;
        revision* target;
    };
    // }}}
    // {{{ revision broadcast get/set
    template<>
    class set<revision, AMBIENT_NUM_PROCS> : public cfunctor {
    public:
        template<class T, int NE> friend class set;
        void* operator new (size_t size, void* placement){ return placement; }
        void* operator new (size_t size){ return ambient::bulk.malloc<sizeof(set)>(); }
        void operator delete (void*, void*){ }
        void operator delete (void* ptr){ }

        static void spawn(revision& r);

        set(revision& r);
        template<int NE> set(set<revision,NE>* s);
        virtual void operator >> (int p);
        virtual bool ready();
        virtual void invoke();
        #ifdef AMBIENT_COMPUTATIONAL_DATAFLOW
        virtual const char* name(){ return "set"; }
        #endif
    private:
        bool evaluated;
        std::vector<bool>* states;
        std::vector<request*>* handles;
        revision* target;
    };
    // }}}
    // {{{ transformable broadcast get/set
    template<>
    class get<transformable, AMBIENT_NUM_PROCS> : public cfunctor {
    public:
        void* operator new (size_t size){ return ambient::bulk.malloc<sizeof(get)>(); }
        void operator delete (void* ptr){ }
        static void spawn(transformable& v);
        get(transformable& v);
        virtual bool ready();
        virtual void invoke();
        #ifdef AMBIENT_COMPUTATIONAL_DATAFLOW
        virtual const char* name(){ return "get"; }
        #endif
    private:
        request* handle;
    };
    template<>
    class set<transformable, AMBIENT_NUM_PROCS> : public cfunctor {
    public:
        void* operator new (size_t size, void* placement){ return placement; }
        void* operator new (size_t size){ return ambient::bulk.malloc<sizeof(set)>(); }
        void operator delete (void*, void*){ }
        void operator delete (void* ptr){ }
        static set& spawn(transformable& v);
        set(transformable& v);
        virtual bool ready();
        virtual void invoke();
        #ifdef AMBIENT_COMPUTATIONAL_DATAFLOW
        virtual const char* name(){ return "set"; }
        #endif
    private:
        bool evaluated;
        std::vector<request*>* handles;
        transformable* target;
    };
    // }}}

} } }

#endif
