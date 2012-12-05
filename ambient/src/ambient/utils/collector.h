#ifndef AMBIENT_UTILS_COLLECTOR
#define AMBIENT_UTILS_COLLECTOR

namespace ambient{

    using ambient::models::velvet::history;

    class collector {
    public:
        struct delete_ptr {
            void operator()( history* element ) const;
            void operator()( void* element ) const;
        };

        collector();
        void push_back(void* o);
        void push_back(history* o);
        void clear();
    private:
        std::vector< history* > str;
        std::vector< void* > raw;
    };
}

#endif

