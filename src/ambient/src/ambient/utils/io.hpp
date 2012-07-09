#ifndef AMBIENT_IO
#define AMBIENT_IO

namespace ambient {

    template<typename T> class future;
    bool verbose();

    class io {
    public:
        std::fstream nullio;
        io() : nullio("/dev/null") { }

        template<class T>
        io& operator<<(future<T> const & obj){
            std::cout << obj.calc_value();
            return *this;
        }

        template<class T>
        io& operator<<(T const & obj){
            std::cout << obj;
            return *this;
        }

        io& operator<<(std::ostream& (*pf)(std::ostream&)){
            std::cout << pf;
            return *this;
        }

        void precision(int p){
            std::cout.precision(p);
        }

        void flush(){
            std::cout.flush();
        }
    };

    class mpio {
    public:
        std::fstream nullio;
        mpio() : nullio("/dev/null") {
            this->v = verbose();
        }

        template<class T>
        mpio& operator<<(future<T> const & obj){
            if(v) std::cout << obj.calc_value();
            else nullio << obj.calc_value();
            return *this;
        }

        template<class T>
        mpio& operator<<(T const & obj){
            if(v) std::cout << obj;
            return *this;
        }

        mpio& operator<<(std::ostream& (*pf)(std::ostream&)){
            if(v) std::cout << pf;
            return *this;
        }

        void precision(int p){
            if(v) std::cout.precision(p);
        }

        void flush(){
            if(v) std::cout.flush();
        }
        bool v;
    };

    extern io cout, cerr;
}

#endif
