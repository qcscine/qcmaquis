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
        class entry {
        public:
            inline entry();
            inline operator char* (){ return (char*)this->data; }
            inline operator double* (){ return (double*)this->data; }
            inline operator std::complex<double>* (){ return (std::complex<double>*)this->data; }
            inline void swap(entry&);
            inline void set_memory(void* memory, size_t bound);
            inline void* get_memory();
            inline bool valid();
            inline bool occupied();
            inline bool trylock();
            inline void unlock();
            inline std::list<cfunctor*>& get_assignments();
            void* header;
            void* data;
            std::list<cfunctor*> assignments;
        };

        inline void embed(void* memory, size_t x, size_t y, size_t bound);
        size_t sid;
        memspec* spec;
        bool clean;
        entry content;
        // layout part //

        inline revision* get_parent(){ return parent; }
        inline revision(memspec*, bool clean = false);
        inline ~revision();
        inline entry& block(size_t x, size_t y);
        inline void add_modifier(sfunctor* m);
        inline std::list<sfunctor*>& get_modifiers();
        inline sfunctor* get_generator();
        inline void set_generator(sfunctor*);
        inline void reset_generator();
        sfunctor* generator;
        std::list<sfunctor*> modifiers;
        revision* parent;
    };

} } }

#endif
