#ifndef AMBIENT_UTILS_COLLECTOR
#define AMBIENT_UTILS_COLLECTOR

#include "ambient/utils/touchstack.h"

namespace ambient{

    using ambient::models::velvet::history;

    class collector {
    public:

        collector(){
            this->size = __cilkrts_get_nworkers();
            this->str = new touchstack< history* >[size];
            this->raw = new touchstack< void* >[size];
        }

       ~collector(){
            delete[] str;
            delete[] raw;
        }

        void push_back(void* o){
            this->raw[__cilkrts_get_worker_number()].push_back(o);
        }

        void push_back(history* o){
            this->str[__cilkrts_get_worker_number()].push_back(o);
        }

        void clear(){
            for(int i = 0; i < size; i++){
                while(!str[i].end_reached())
                    delete str[i].pick();
                while(!raw[i].end_reached())
                    ambient::static_memory::free<FUTURE_SIZE>(raw[i].pick());
            }
        }

    private:
        touchstack< history* >* str;
        touchstack< void* >* raw;
        int size;
    };

}

#endif
