#ifndef AMBIENT_INTERFACE_FUTURE
#define AMBIENT_INTERFACE_FUTURE
// see history for an advanced version // supports multiple revisions

namespace ambient {

    template <typename T>
    class future {
    private:
        future(){}
        template<typename S> inline future& operator = (const S& v){ }
        template<typename S> inline future& operator = (const future& v){ }
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
           *value = (T)f; // move semantics
        }

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

        inline operator T () const {
            if(value == NULL){
                ambient::playout();
                value = (T*)&(*ghost);
            }
            return *value;
        }

        inline const T& get_value() const {
            return *(T*)&(*ghost);
        }

        inline T& get_value(){
            return *(T*)&(*ghost);
        }

        inline future<T>& unfold(){ // should be called reset
            this->value = NULL;
            return *this;
        }

        inline const future<T>& unfold() const {
            return *this;
        }
        ptr    ghost;
    private:
        mutable T* value;
    };

    template<typename T1, typename T2>
    inline const T2 operator / (T1 lhs, future<T2> rhs){ 
        return (lhs / (T2)rhs); 
    }

    template<typename T>
    inline const future<T> operator + (future<T> lhs, future<T> rhs){ 
        return future<T>((T)lhs + (T)rhs); // explicit
    }

}

#endif
