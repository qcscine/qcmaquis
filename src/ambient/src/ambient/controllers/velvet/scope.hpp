#ifndef AMBIENT_INTERFACE_SCOPE
#define AMBIENT_INTERFACE_SCOPE

namespace ambient { 

    using ambient::controllers::velvet::controller;

    template<scope_t T>
    class scope {};

    template<>
    class scope<base> : public controller::scope {
    public:
        scope(){
            this->tunable = true;
            this->round = ambient::channel.volume;
            this->state = ambient::rank() ? REMOTE : LOCAL;
            this->sector = 0;
        }
        /*class info {
        public:
            info(){
                footprint = remote = local =
                pin = load[0] = load[1] = rank = 0;
            }
            void repeat(){
                pin = remote = local =
                footprint = 0;
            }
            void clear(){
                load[0] = load[1] = 0;
            }
            int decide(){
                if(power != complexity::N3 && local != remote){
                    if(local > remote){ load[ambient::rank()] += footprint; return ambient::rank(); }
                    else{ load[1-ambient::rank()] += footprint; return (1-ambient::rank()); }
                }

                if(load[0] / std::max(1,load[1]) > 2) rank = 1;
                if(load[1] / std::max(1,load[0]) > 2) rank = 0;
                load[rank] += footprint; //std::pow(pin, power);

                return rank;
            }
            void add_as_new(size_t size){
                footprint += size;
                pin = std::max(pin, size);
            }
            void add_as_local(size_t size){
                pin = std::max(pin, size);
                local += size;
            }
            void add_as_remote(size_t size){
                pin = std::max(pin, size);
                remote += size;
            }
            size_t footprint;
            size_t remote;
            size_t local;
            size_t pin;
            size_t power;
            int load[2];
            int rank;
        } tuning;*/
    };

    template<>
    class scope<single> : public controller::scope {
    public:
        static int iterator;
        static int factor;
        static int effect;
        static void compact(size_t n){
            iterator = 0;
            if(n <= ambient::channel.volume) return;
            factor = (int)(n / ambient::channel.volume); // iterations before switch
            effect = (int)n;
        }
        scope(){
            this->tunable = false;
            ambient::controller.set_context(this);
            this->round = ambient::channel.volume;
            this->sector = (++iterator %= round*factor)/factor;
            this->state = (this->sector == ambient::rank()) ? LOCAL : REMOTE;
        }
       ~scope(){
            if(effect && !--effect) factor = 1;
            ambient::controller.pop_context();
        }
    };

    template<>
    class scope<shared> : public controller::scope {
    public:
        scope(){
            this->tunable = false;
            ambient::controller.set_context(this);
            this->round = ambient::channel.volume;
            this->state = COMMON;
        }
       ~scope(){
            ambient::controller.pop_context();
        }
    };
}

#endif
