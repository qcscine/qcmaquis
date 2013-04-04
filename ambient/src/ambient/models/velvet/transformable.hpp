namespace ambient { namespace models { namespace velvet {

    template<typename T> constexpr T op_single (T a){ }
    template<typename T> constexpr T op_double (T a, T b){ }

    template<typename T> constexpr T op_sqrt (T a) { return std::sqrt(a); }
    template<typename T> constexpr T op_plus (T a, T b){ return a + b; }
    template<typename T> constexpr T op_minus(T a, T b){ return a - b; }
    template<typename T> constexpr T op_mul  (T a, T b){ return a * b; }
    template<typename T> constexpr T op_div  (T a, T b){ return a / b; }

    inline transformable::numeric_union::operator bool& (){ return b; }
    inline transformable::numeric_union::operator double& (){ return d; }
    inline transformable::numeric_union::operator std::complex<double>& (){ return c; }
    inline void transformable::numeric_union::operator = (bool value) { b = value; }
    inline void transformable::numeric_union::operator = (double value) { d = value; }
    inline void transformable::numeric_union::operator = (std::complex<double> value) { c = value; }

    inline void* transformable::operator new (size_t size, void* placement){ return placement; }
    inline void transformable::operator delete (void*, void*){ /* doesn't throw */ }

    template<typename T>
    inline transformable_value<T>::transformable_value(T value){
        ambient::model.index(this);
        this->v = value;
    }

    template<typename T>
    inline transformable::numeric_union transformable_value<T>::eval() const { 
        return this->v; 
    }

    template<typename T>
    inline transformable& transformable_value<T>::operator += (transformable& r){ 
        return *new(this)transformable_expr<T, decltype(op_plus<T>), op_plus>(
                   &r, 
                   new (ambient::pool.malloc<AMBIENT_FUTURE_SIZE>()) transformable_value(*this)
               );
    }

    template <typename T, typename FP, FP OP>
    inline transformable::numeric_union transformable_expr<T,FP,OP>::eval() const {
        this->v = OP(this->l->eval(), this->r->eval());
        ambient::pool.free<AMBIENT_FUTURE_SIZE>((void*)this->l);
        ambient::pool.free<AMBIENT_FUTURE_SIZE>((void*)this->r);
        new((void*)this)transformable_value<T>(*(transformable_value<T>*)this);
        ambient::model.index((transformable*)this);
        return this->v;
    }

    template <typename T, typename FP, FP OP>
    inline transformable& transformable_expr<T,FP,OP>::operator += (transformable& r){ 
        return *new(this)transformable_expr<T, decltype(op_plus<T>), op_plus>(
                   &r, 
                   new (ambient::pool.malloc<AMBIENT_FUTURE_SIZE>()) transformable_expr<T,FP,OP>(*this)
               ); 
    }

    template <typename T, typename FP, FP OP>
    inline transformable_expr<T,FP,OP>::transformable_expr(const transformable* l){ 
        this->l = l; 
    }

    template <typename T, typename FP, FP OP>
    inline transformable_expr<T,FP,OP>::transformable_expr(const transformable* l, const transformable* r){ 
        this->l = l; 
        this->r = r; 
    }

} } }

