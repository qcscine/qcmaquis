#ifndef AMBIENT_UTILS_TIMINGS
#define AMBIENT_UTILS_TIMINGS
#include "ambient/ambient.hpp"
#include "ambient/utils/io.hpp"
#include <pthread.h>

#define BILLION 0x3B9ACA00

namespace ambient {

    class timer {
    public:
        timer(std::string name, pthread_t thread): val(0.0), name(name), count(0), thread_(thread){}
        timer(std::string name): val(0.0), name(name), count(0), thread_(pthread_self()){}
       ~timer(){ report(); }
     
        double get_time() const {
            return val;
        }
        void report(){
            ambient::cout << name << " " << val << ", count : " << count << "\n";
        }
        void begin(){
            pthread_getcpuclockid(this->thread_,&this->cid_);
            struct timespec ts; // from time.h
            clock_gettime(this->cid_, &ts);
            this->t0 = ts.tv_sec+(((double)ts.tv_nsec / (double)(BILLION)));
        }    
        void end(){
            count++;
            struct timespec ts;
            clock_gettime(this->cid_, &ts);
            this->val += (ts.tv_sec+(((double)ts.tv_nsec / (double)(BILLION))) - this->t0);
        }
    private:    
        pthread_t thread_; 
        clockid_t cid_;
        double val, t0;
        unsigned long long count;
        std::string name;
    };
}

#endif

