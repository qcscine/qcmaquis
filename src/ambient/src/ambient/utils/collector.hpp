#ifndef AMBIENT_UTILS_COLLECTOR
#define AMBIENT_UTILS_COLLECTOR

#define COLLECTOR_STR_RESERVE 65536
#define COLLECTOR_RAW_RESERVE 1024

namespace ambient{

    using ambient::models::velvet::history;

    class collector {
    public:

        collector(){
            this->str.reserve(COLLECTOR_STR_RESERVE);
            this->raw.reserve(COLLECTOR_RAW_RESERVE);
        }

        void push_back(void* o){
            this->raw.push_back(o);
        }

        void push_back(history* o){
            this->str.push_back(o);
        }

        struct delete_ptr {
            void operator()( history* element ) const {
                delete element;
            }
            void operator()( void* element ) const {
                ambient::static_memory::free<FUTURE_SIZE>(element);
            } 
        };

        void clear(){
            std::for_each( str.begin(), str.end(), delete_ptr());
            std::for_each( raw.begin(), raw.end(), delete_ptr());
            str.clear();
            raw.clear();
        }

    private:
        std::vector< history* > str;
        std::vector< void* > raw;
    };

}

#endif
