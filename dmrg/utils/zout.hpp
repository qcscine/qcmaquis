#ifndef DMRG_COUT_H
#define DMRG_COUT_H

#ifdef MPI_PARALLEL
#include "ambient/ambient.h"
#endif
#include <iostream>
#include <fstream>

class master_cout {
public:

    std::fstream nullio;
    master_cout():nullio("/dev/null"){}

    template<class T>
    master_cout& operator <<(T const & obj)
    {
#ifdef MPI_PARALLEL
        !ambient::is_master() ? this->nullio << obj :
#endif
        std::cout << obj;
        return *this;
    }
    master_cout& operator <<(std::ostream& (*pf)(std::ostream&))
    {
#ifdef MPI_PARALLEL
        !ambient::is_master() ? this->nullio << pf :
#endif
        std::cout << pf;
        return *this;
    }
    void precision(int p)
    {
        std::cout.precision(p);
    }
};

class group_master_cout : public master_cout
{
public:
    template<class T>
    group_master_cout& operator <<(T const & obj)
    {
#ifdef MPI_PARALLEL
        !ambient::is_group_master() ? this->nullio << obj :
#endif
        std::cout << obj;
        return *this;
    }

    group_master_cout& operator <<(std::ostream& (*pf)(std::ostream&))
    {
#ifdef MPI_PARALLEL
        !ambient::is_group_master() ? this->nullio << pf :
#endif
        std::cout << pf;
        return *this;
    }
};

extern master_cout zout;
extern group_master_cout gzout;

#endif
