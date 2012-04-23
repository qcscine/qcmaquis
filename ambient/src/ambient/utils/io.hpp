#ifndef AMBIENT_IO_HPP
#define AMBIENT_IO_HPP

#include <iostream>
#include <fstream>

namespace ambient {

    bool verbose(){
        return (rank()?false:true); 
    }

    class io {
    public:
        std::fstream nullio;
        io()
        : nullio("/dev/null")
        {
        }

        template<class T>
        io& operator<<(T const & obj){
            if(verbose()) std::cout << obj;
            return *this;
        }

        io& operator<<(std::ostream& (*pf)(std::ostream&)){
            if(verbose()) std::cout << pf;
            return *this;
        }

        void precision(int p){
            if(verbose()) std::cout.precision(p);
        }
    } cout;
}

#endif
