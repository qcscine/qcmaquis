#ifndef AMBIENT_INTERFACE_FUTURE
#define AMBIENT_INTERFACE_FUTURE
// see history for an advanced version // supports multiple revisions

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
#endif

//#define FUTURE_SAFE_CHECK

#ifndef RVALUE
#define RVALUE
#endif

namespace ambient {

    template <typename T>
    class future {
    private:
        future(){}
        template<typename S> inline future& operator = (const S& v){ }
    public:
        typedef typename boost::intrusive_ptr< container<sizeof(T)> > ptr;
        typedef T value_type;

        explicit inline future(const ptr& p)
        : ghost(p), value((T*)&(*p))
        {
        }

        explicit inline future(const future& f){
            ghost = new container<sizeof(T)>();
            value = (T*)&(*ghost);
           *value = f.calc_value();
        }

        template<typename S>
        inline future(const future<S>& f){ // can be optimized later
            ghost = new container<sizeof(T)>();
            value = (T*)&(*ghost);
           *value = (T)f.calc_value();
        }

        inline future& operator = (const future& v){ 
            ghost = v.ghost;
            value = (T*)&(*ghost);
            return *this;
        }

#ifdef RVALUE
        inline future(future&& f){
            printf("RVALUE COPY!\n");
            ghost = f.ghost;
            value = (T*)&(*ghost);
        }

        inline future& operator = (future&& v){ 
            printf("RVALUE =!\n");
            ghost = v.ghost;
            value = (T*)&(*ghost);
            return *this;
        }
#endif

        inline future(double v){
            ghost = new container<sizeof(T)>();
            value = (T*)&(*ghost);
           *value = v;
        }

        inline future(std::complex<double> v){
            ghost = new container<sizeof(T)>();
            value = (T*)&(*ghost);
           *value = v;
        }

        inline T calc_value() const {
            if(value == NULL){
                ambient::playout();
                value = (T*)&(*ghost);
            }
            return *value;
        }

        inline operator T () const {
            return this->calc_value();
        }

        inline const T& get_value() const {
            return *(T*)&(*ghost);
        }

        inline T& get_value(){
            return *(T*)&(*ghost);
        }

        inline future<T>& unfold(){ // should be called reset
#ifdef FUTURE_SAFE_CHECK
            if(this->value == NULL) printf("ERROR: Overusing the future!\n");
#endif
            this->value = NULL;
            return *this;
        }

        inline const future<T>& unfold() const {
            return *this;
        }

#ifdef HAVE_ALPS_HDF5
        inline void load(alps::hdf5::archive & ar) { /*ambient::cerr << "I don't do much." << std::endl;*/ }
        inline void save(alps::hdf5::archive & ar) const { /*ambient::cerr << "I don't do much either." << std::endl;*/ }
#endif

        ptr    ghost;
    private:
        mutable T* value;
    };

    template<typename T>
    inline const T operator / (double lhs, const future<T>& rhs){ 
        return (lhs / rhs.calc_value()); 
    }

    template<typename T>
    inline const T operator / (std::complex<double> lhs, const future<T>& rhs){ 
        return (lhs / rhs.calc_value()); 
    }

    template<typename T>
    inline T operator / (const future<T>& lhs, const future<T>& rhs){ 
        return (lhs.calc_value() / rhs.calc_value()); 
    }

    template<typename T>
    inline const future<T> operator + (const future<T>& lhs, const future<T>& rhs){ 
        return future<T>(lhs.calc_value() + rhs.calc_value()); // explicit
    }

    inline double sqrt(const future<double>& f){
        return std::sqrt(f.calc_value());
    } 

}

namespace alps { namespace numeric {

    inline double real(double f){
        return f;
    }

    inline double real(std::complex<double> f){
        return f.real();
    }

    inline const ambient::future<double>& real(const ambient::future<double>& f){
        return f;
    }

    inline double real(const ambient::future<std::complex<double> >& f){
        return f.calc_value().real();
    }

    inline const std::vector<double>& real(const std::vector<double>& f){
        return f;
    }

    inline std::vector<double> real(const std::vector<std::complex<double> >& f){
        std::vector<double> res;
        for(size_t k = 0; k < f.size(); ++k) res.push_back(f[k].real());
        return res;
    }

    inline std::vector<double> real(const std::vector<ambient::future<double> >& f){
        ambient::playout();
        std::vector<double> res;
        for(size_t k = 0; k < f.size(); ++k) res.push_back(f[k].get_value());
        return res;
    }

    inline std::vector<double> real(const std::vector<ambient::future<std::complex<double> > >& f){
        ambient::playout();
        std::vector<double> res;
        for(size_t k = 0; k < f.size(); ++k) res.push_back(f[k].get_value().real());
        return res;
    }

} }

#endif
