#ifndef AMBIENT_ZOUT_H
#define AMBIENT_ZOUT_H

#ifdef MPI_PARALLEL
#include "ambient/ambient.h"
#endif
#include <iostream>
#include <fstream>

namespace ambient {

    bool outlet();
    bool goutlet(); 
    // @ambient.hpp

    class master_cout {
    public:
        std::fstream nullio;
        master_cout()
        : nullio("/dev/null")
        {
        }

        template<class T>
        master_cout& operator <<(T const & obj){
#ifdef MPI_PARALLEL
            !ambient::outlet() ? this->nullio << obj :
#endif
            std::cout << obj;
            return *this;
        }

        master_cout& operator <<(std::ostream& (*pf)(std::ostream&)){
#ifdef MPI_PARALLEL
            !ambient::outlet() ? this->nullio << pf :
#endif
            std::cout << pf;
            return *this;
        }

        void precision(int p){
            std::cout.precision(p);
        }
    };

    class group_master_cout : public master_cout
    {
    public:
        template<class T>
        group_master_cout& operator <<(T const & obj){
#ifdef MPI_PARALLEL
            !ambient::goutlet() ? this->nullio << obj :
#endif
            std::cout << obj;
            return *this;
        }

        group_master_cout& operator <<(std::ostream& (*pf)(std::ostream&)){
#ifdef MPI_PARALLEL
            !ambient::goutlet() ? this->nullio << pf :
#endif
            std::cout << pf;
            return *this;
        }
    };

    extern master_cout zout;
    extern group_master_cout gzout;

}

#endif
