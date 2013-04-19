#ifndef MAQUIS_TIMINGS_H
#define MAQUIS_TIMINGS_H

#include <string>
#include <fstream>
#include <iostream>
#include <boost/chrono.hpp>
#include "utils/io.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif

class Timer
{
public:
    Timer(std::string name_)
    : val(0.0), name(name_), nCounter(0) { }
    
    ~Timer() { maquis::cout << name << " " << val << ", nCounter : " << nCounter << std::endl; }
    
    Timer & operator+=(double t) {
        val += t;
        return *this;
    }
    
    void begin() {
        t0 = boost::chrono::system_clock::now();
    }
    
    void end() {
		nCounter += 1;
        boost::chrono::duration<double> sec = boost::chrono::system_clock::now() - t0;
        val += sec.count();
    }
    
    double get_time() const {
	    return  val;
    }
    
    friend std::ostream& operator<< (std::ostream& os, Timer const& timer) {
        os << timer.name << " " << timer.val << ", nCounter : " << timer.nCounter;
        return os;
    }
    
protected:
    double val;
    boost::chrono::system_clock::time_point t0;
    unsigned long long nCounter;
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
