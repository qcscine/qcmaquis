#ifndef AMBIENT_UTILS_TIMINGS
#define AMBIENT_UTILS_TIMINGS

#define BILLION 0x3B9ACA00

namespace ambient {

    class timer {
    public:
        timer(std::string name, pthread_t thread): val(0.0),name(name),count(0),thread_(thread){}
        timer(std::string name): val(0.0),name(name),count(0),thread_(pthread_self()){}
        ~timer(){ std::cout << name << " " << val << ", count : " << count << std::endl; }
    
        inline void begin(){
             pthread_getcpuclockid(this->thread_,&this->cid_);
             struct timespec ts; //from time.h
             clock_gettime(this->cid_, &ts);
             this->t0 = ts.tv_sec+(((double)ts.tv_nsec / (double)(BILLION)));
        }    
        
        inline void end(){
            count ++;
            struct timespec ts; //from time.h
            clock_gettime(this->cid_, &ts);
            this->val += (ts.tv_sec+(((double)ts.tv_nsec / (double)(BILLION))) - this->t0);
        }
    
        inline double get_time() const { return val; }    
       
        inline timer & operator+=(double t) {
            val += t;
            return *this;
        }
     
        friend std::ostream& operator<< (std::ostream& os, timer const& timer) {
            os << timer.name << " " << timer.val << ", count : " << timer.count;
            return os;
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

