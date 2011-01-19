#ifndef TIMINGS_H
#define TIMINGS_H

unsigned long long getcpuclocks() {
    unsigned long long clk;
    __asm { rdtsc };
    __asm { shl rdx,32 };
    __asm { add rax,rdx };
    __asm { mov clk,rax };
    return clk;
}

#define CPU_FREQ 2.2e9

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
