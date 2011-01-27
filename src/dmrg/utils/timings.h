#ifndef TIMINGS_H
#define TIMINGS_H

#include "utils/zout.hpp"

#ifdef __CUBLAS__
#include "cuda_runtime_api.h" 
#endif


/*
unsigned long long getcpuclocks() {
    unsigned long long clk;
    __asm__ ("rdtsc");
    __asm__("shl %rdx,32");
    __asm { add rax,rdx };
    __asm { mov clk,rax };
    return clk;
}*/


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
    
    ~Timer() { zout << name << " " << val << ", nCounter : " << nCounter << std::endl; }
    
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
            1+1;
        else
            val += (getcpuclocks()-t0)/freq;
    }

    const double GetTime()
    {
	return  val;
    }    
protected:
    double val, t0;
    unsigned long long freq, nCounter;
    std::string name;
};

#ifdef __CUBLAS__

class TimerCuda : public Timer
{
public:
	TimerCuda(std::string name_) : Timer(name_)
	{}
	
	~TimerCuda(){}
	
	void begin()
    {
        cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord(start, 0);
    }
    
    void end()
    {
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		float elapsedTime;
		cudaEventElapsedTime(&elapsedTime, start, stop); // that's our time!
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
		val = static_cast<double>(elapsedTime)/1000 ; // because time is in microsecond
    }
	
	
private:
	cudaEvent_t start, stop;
	
};

#endif


#endif
