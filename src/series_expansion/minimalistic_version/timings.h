#ifndef TIMINGS_H
#define TIMINGS_H
#include "vli/vli_config.h"
#include <string>
#include <fstream>
#ifdef _OPENMP
#include <omp.h>
#endif //_OPENMP


#ifdef __CUBLAS__
#include "cuda_runtime_api.h"
#endif

//prototype for removing warning of my compiler
unsigned long getcpuclocks();

#if defined(__i386__)
 
unsigned long getcpuclocks()
{
 unsigned long tsc;
 asm(".byte 0x0f, 0x31" : "=A" (x));
 return( tsc );
}
 
#elif (defined(__amd64__) || defined(__x86_64__))
 
unsigned long getcpuclocks()
{
 unsigned long lo, hi;
 asm( "rdtsc" : "=a" (lo), "=d" (hi) );
 return( lo | (hi << 32) );
}

#endif


class Timer
{
public:
    Timer(std::string name_)
    : val(0), name(name_), freq(CPU_FREQ),nCounter(0)
    { }
    
    ~Timer() { std::cout << name << " " << val << ", nCounter : " << nCounter << std::endl; }
    
    Timer & operator+=(double t)
    {
        val += t;
        return *this;
    }
    
    void begin()
    {
        t0 = getcpuclocks();
    }
    
    void end()
    {
		nCounter += 1;
        unsigned long long t1 = getcpuclocks();
        if (t1 < t0)
            assert(true);
        else
            val += (getcpuclocks()-t0)/freq; // to milliseconds
    }

    const double GetTime()
    {
	return  val;
    }    

    void save()
    {
       std::ofstream o;
       std::string file("time");
       file+=name; 
       o.open(file.c_str(),std::ios::app);
#ifdef _OPENMP 
           o << val << " "  <<  omp_get_max_threads() <<  std::endl;
#else
           o << val << " "  <<  1 <<   std::endl;
#endif
       o.close();
    }
protected:
    double val, t0;
    std::string name;
    unsigned long long freq, nCounter;
};

#endif
