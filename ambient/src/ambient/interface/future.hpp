#ifndef AMBIENT_INTERFACE_FUTURE
#define AMBIENT_INTERFACE_FUTURE

#include <alps/numeric/real.hpp> 

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
        template<typename S> inline future& operator = (const S& v){ }
    public:
        typedef void* ptr;
        typedef T value_type;

        future()
        : symlink(false)
        {
            printf("ambient::future$ used default constructor\n");
            ghost = ambient::pool.malloc<FUTURE_SIZE>();
            value = (T*)ghost;
           *value = T();
        }

        explicit inline future(const ptr& p)
        : ghost(p), value((T*)p), symlink(true)
        {
        }

        inline ~future(){
            if(!symlink) ambient::destroy(ghost);
        }

        explicit inline future(const future& f)
        : symlink(false) 
        {
            ghost = ambient::pool.malloc<FUTURE_SIZE>();
            value = (T*)ghost;
           *value = f.calc_value();
        }

        template<typename S>
        inline future(const future<S>& f) // can be optimized later
        : symlink(false)
        {
            const_cast<future<S>&>(f).symlink = true;
            ghost = f.ghost;
            value = (T*)ghost;
        }

        /*template<typename S>
        inline future(const future<S>& f) // can be optimized later
        : symlink(false) 
        {
            ghost = ambient::pool.malloc<FUTURE_SIZE>();
            value = (T*)ghost;
           *value = (T)f.calc_value();
        }*/

        inline future& operator = (const future& f){
            const_cast<future&>(f).symlink = true;
            ghost = f.ghost;
            value = (T*)ghost;
            return *this;
        }

#ifdef RVALUE
        inline future(future&& f){
            f.symlink = true;
            ghost = f.ghost;
            value = (T*)ghost;
        }

        inline future& operator = (future&& f){
            f.symlink = true;
            ghost = f.ghost;
            value = (T*)ghost;
            return *this;
        }
#endif

        inline future(double v)
        : symlink(false)
        {
            ghost = ambient::pool.malloc<FUTURE_SIZE>();
            value = (T*)ghost;
           *value = v;
        }

        inline future(std::complex<double> v)
        : symlink(false)
        {
            ghost = ambient::pool.malloc<FUTURE_SIZE>();
            value = (T*)ghost;
           *value = v;
        }

        inline T calc_value() const {
            if(value == NULL){
                ambient::sync();
                value = (T*)ghost;
            }
            return *value;
        }

        inline future& operator += (const future& rhs){ 
            *(T*)ghost += rhs.calc_value();
            return *this;
        }

        inline operator T () const {
            return this->calc_value();
        }

        inline const T& get_value() const {
            return *(T*)ghost;
        }

        inline T& get_value(){
            return *(T*)ghost;
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
        bool   symlink;
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

    template<typename T>
    inline const ambient::future<double>& real(const ambient::future<T>& f){
        return *(ambient::future<double>*)&f;
    }

    template<typename T>
    inline std::vector<double> real(const std::vector<ambient::future<T> >& f){
        ambient::sync();
        std::vector<double> res;
        for(size_t k = 0; k < f.size(); ++k) res.push_back(*(double*)&f[k].get_value());
        return res;
    }

}

#endif
