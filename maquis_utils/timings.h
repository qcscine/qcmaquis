#ifndef MAQUIS_TIMINGS_H
#define MAQUIS_TIMINGS_H

#include <string>
#include <fstream>
#include <iostream>
#include "utils/io.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif

#ifdef __CUBLAS__
#include "cuda_runtime_api.h" 
#endif

#define BILLION 0x3B9ACA00

#if defined(__i386__)
 
inline unsigned long getcpuclocks()
{
 unsigned long tsc;
 asm(".byte 0x0f, 0x31" : "=A" (tsc));
 return( tsc );
}
 
#elif (defined(__amd64__) || defined(__x86_64__))
 
inline unsigned long getcpuclocks()
{
 unsigned long lo, hi;
 asm( "rdtsc" : "=a" (lo), "=d" (hi) );
 return( lo | (hi << 32) );
}

#elif (defined(__powerpc64__) || defined(__powerpc__))

inline unsigned long getcpuclocks()
{
  unsigned long long int result=0;
  unsigned long int upper, lower,tmp;
  __asm__ volatile(
                "0:                  \n"
                "\tmftbu   %0           \n"
                "\tmftb    %1           \n"
                "\tmftbu   %2           \n"
                "\tcmpw    %2,%0        \n"
                "\tbne     0b         \n"
                : "=r"(upper),"=r"(lower),"=r"(tmp)
                );
  result = upper;
  result = result<<32;
  result = result|lower;
  return(result);
}

#else

inline unsigned long getcpuclocks() { return 0; }

#endif

class Timer
{
public:
    Timer(std::string name_)
    : val(0.0), name(name_), freq((long long unsigned int)CPU_FREQ), nCounter(0) { }
    
    ~Timer() { maquis::cout << name << " " << val << ", nCounter : " << nCounter << std::endl; }
    
    Timer & operator+=(double t) {
        val += t;
        return *this;
    }
    
    void begin() {
        t0 = getcpuclocks();
    }
    
    void end() {
		nCounter += 1;
        unsigned long long t1 = getcpuclocks();
        if (t1 > t0)
            val += (getcpuclocks()-t0)/freq;
    }

    double get_time() const {
	    return  val;
    }    
  
    friend std::ostream& operator<< (std::ostream& os, Timer const& timer) {
        os << timer.name << " " << timer.val << ", nCounter : " << timer.nCounter;
        return os;
    }
    
protected:
    double val, t0;
    unsigned long long freq, nCounter;
    std::string name;
};

#ifdef _OPENMP
class TimerOMP : public Timer {
public:
	TimerOMP(std::string name_) : Timer(name_), timer_start(0.0), timer_end(0.0){}

	~TimerOMP(){}
	
	void begin() {
		timer_start = omp_get_wtime(); 
	}
	
	void end() {
		timer_end = omp_get_wtime();
		val += timer_end - timer_start;
	}
private:
	double timer_start, timer_end;
};
#endif

#ifdef AMBIENT 
class TimerPTH : public Timer{
public:
    TimerPTH(std::string name, pthread_t thread): Timer(name),thread_(thread){}
    TimerPTH(std::string name): Timer(name),thread_(pthread_self()){}
    ~TimerPTH(){}

    void begin(){
         pthread_getcpuclockid(thread_,&cid_);
    }    
  
    void end(){
        nCounter ++;
        struct timespec ts; //from time.h
        clock_gettime(cid_, &ts);
        val += ts.tv_sec+(((double)ts.tv_nsec / (double)(BILLION)));
    }
private:    
    pthread_t thread_; 
    clockid_t cid_;
};
#endif

#ifndef WIN32
#include <sys/time.h>
#else

#include <ctime>

struct timeval {
	time_t tv_sec, tv_usec;
};

void gettimeofday(timeval * tv, void *)
{
	tv->tv_sec = time(NULL);
	tv->tv_usec = 0;
}

#endif

#endif
