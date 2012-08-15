#ifndef AMBIENT_UTILS_COLLECTOR
#define AMBIENT_UTILS_COLLECTOR

#include "ambient/utils/touchstack.h"

namespace ambient{

    using ambient::models::velvet::history;

    class collector {
    public:

        collector(){
            this->str = new touchstack< history* >[__cilkrts_get_nworkers()];
            this->raw = new touchstack< void* >[__cilkrts_get_nworkers()];
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
            for(int i = 0; i < __cilkrts_get_nworkers(); i++){
                str[i].clear();
                raw[i].purge();
            }
        }

    private:
        touchstack< history* >* str;
        touchstack< void* >* raw;
    };

}

#endif
