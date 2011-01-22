#ifndef DMRG_COUT_H
#define DMRG_COUT_H

#ifdef MPI_PARALLEL
#include "ambient/ambient.h"
#endif

class dmrg_cout {
public:
    template<class T>
    dmrg_cout& operator <<(T const & obj)
    {
#ifdef MPI_PARALLEL
        if(ambient::scheduler::instance().is_ambient_master())
#endif
        std::cout << obj;
        return *this;
    }

    dmrg_cout& operator <<(std::ostream& (*pf)(std::ostream&))
    {
#ifdef MPI_PARALLEL
        if(ambient::scheduler::instance().is_ambient_master())
#endif
        std::cout << pf;
        return *this;
    }

    void precision(int p)
    {
        std::cout.precision(p);
    }
} zout;

#endif
