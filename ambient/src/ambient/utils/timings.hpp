#ifndef AMBIENT_UTILS_TIMINGS
#define AMBIENT_UTILS_TIMINGS

#define BILLION 0x3B9ACA00

class __a_timer {
public:
    __a_timer(std::string name, pthread_t thread): val(0.0),name(name),nCounter(0),thread_(thread){}
    __a_timer(std::string name): val(0.0),name(name),nCounter(0),thread_(pthread_self()){}
    ~__a_timer(){ std::cout << name << " " << val << ", nCounter : " << nCounter << std::endl; }

    inline void begin(){
         pthread_getcpuclockid(this->thread_,&this->cid_);
         struct timespec ts; //from time.h
         clock_gettime(this->cid_, &ts);
         this->t0 = ts.tv_sec+(((double)ts.tv_nsec / (double)(BILLION)));
    }    
    
    inline void end(){
        nCounter ++;
        struct timespec ts; //from time.h
        clock_gettime(this->cid_, &ts);
        this->val += (ts.tv_sec+(((double)ts.tv_nsec / (double)(BILLION))) - this->t0);
    }

    inline double get_time() const { return val; }    
   
    inline __a_timer & operator+=(double t) {
        val += t;
        return *this;
    }
 
    friend std::ostream& operator<< (std::ostream& os, __a_timer const& timer) {
        os << timer.name << " " << timer.val << ", nCounter : " << timer.nCounter;
        return os;
    }

private:    
    pthread_t thread_; 
    clockid_t cid_;
    double val, t0;
    unsigned long long nCounter;
    std::string name;
};

#endif

