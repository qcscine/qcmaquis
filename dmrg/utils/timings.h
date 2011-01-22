#ifndef TIMINGS_H
#define TIMINGS_H

/*
unsigned long long getcpuclocks() {
    unsigned long long clk;
    __asm__ ("rdtsc");
    __asm__("shl %rdx,32");
    __asm { add rax,rdx };
    __asm { mov clk,rax };
    return clk;
}*/

#define CPU_FREQ 1e9

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
    : val(0), name(name_), freq(CPU_FREQ)
    { }
    
    ~Timer() { cout << name << " " << val << endl; }
    
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
        unsigned long long t1 = getcpuclocks();
        if (t1 < t0)
            1+1;
        else
            val += (getcpuclocks()-t0)/freq;
    }
    
private:
    double val, t0;
    unsigned long long freq;
    std::string name;
};

#endif
