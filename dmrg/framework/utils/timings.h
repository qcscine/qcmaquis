#ifndef MAQUIS_TIMINGS_H
#define MAQUIS_TIMINGS_H

#include <string>
#include <fstream>
#include <iostream>
#include <chrono>
#include <utility>
#include "utils/io.hpp"

#ifdef MAQUIS_OPENMP
#include "omp.h"
#endif

class Timer
{
public:
    Timer(std::string  name_) : name(std::move(name_)) {}

    ~Timer() { maquis::cout << name << " took " << val << " [s], nCounter : " << nCounter << std::endl; }

    Timer & operator+=(double t) {
        val += t;
        return *this;
    }

    virtual void begin() {
        t0 = std::chrono::high_resolution_clock::now();
    }

    void end() {
      nCounter += 1;
      std::chrono::duration<double> sec = std::chrono::high_resolution_clock::now() - t0;
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
    double val = 0.0;
    std::string name;
    std::chrono::time_point<std::chrono::high_resolution_clock> t0;
    unsigned long long nCounter = 0;
};

#ifdef MAQUIS_OPENMP
class TimerOMP : public Timer {
public:
    TimerOMP(std::string name_) : Timer(name_) {}

    ~TimerOMP() = default;

    void begin() {
        timer_start = omp_get_wtime();
    }

    void end() {
        timer_end = omp_get_wtime();
        val += timer_end - timer_start;
    }
private:
    double timer_start = 0.0;
    double timer_end = 0.0;
};
#endif


#endif
