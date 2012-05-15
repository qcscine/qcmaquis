#ifndef AMBIENT_INTERFACE_FUTURE_T
#define AMBIENT_INTERFACE_FUTURE_T
// see history for an advanced version // supports multiple revisions

namespace ambient {

    template <typename T>
    class future {
    public:
        typedef typename boost::intrusive_ptr< container<sizeof(T)> > ptr;
        typedef T value_type;

        future()
        : value(NULL)
        {
            this->naked = new container<sizeof(T)>();
            this->ghost = (container<sizeof(T)>*)this->naked;
        }

        future(const future& f){
            this->naked = new container<sizeof(T)>();
            *(T*)this->naked = (T)f; // unfolding f (can be more optimal / wo playout)
            this->ghost = (container<sizeof(T)>*)this->naked;
            this->value = NULL; // let's playout again for copy
        }

        future(double value){
            this->naked = new container<sizeof(T)>();
            this->ghost = (container<sizeof(T)>*)this->naked;
            *(T*)this->naked = value;
            this->value = NULL; //(T*)this->naked; // let's playout
        }

        future(std::complex<double> value){
            this->naked = new container<sizeof(T)>();
            this->ghost = (container<sizeof(T)>*)this->naked;
            *(T*)this->naked = value;
            this->value = NULL; //(T*)this->naked; // let's playout
        }

        operator T () const {
            if(this->value == NULL){
                ambient::playout();
                this->value = (T*)&(*this->ghost);
            }
            return *this->value;
        }

        inline T*& unfold(){
            return (T*&)this->naked;
        }

        inline const T*& unfold() const {
            return (const T*&)this->naked;
        }
    private:
        ptr    ghost;
        mutable T* value;
        void*  naked;
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
