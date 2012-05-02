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

        future(double value)
        {
            this->naked = new container<sizeof(T)>();
            this->ghost = (container<sizeof(T)>*)this->naked;
            this->value = (T*)this->naked;
            *this->value = value;
        }

        future(std::complex<double> value)
        {
            this->naked = new container<sizeof(T)>();
            this->ghost = (container<sizeof(T)>*)this->naked;
            this->value = (T*)this->naked;
            *this->value = value;
        }

        operator T (){
            if(this->value == NULL){
                ambient::playout();
                this->value = (T*)&(*this->ghost);
            }
            return *this->value;
        }

        T*& unfold() {
            return (T*&)this->naked;
        }

        const T*& unfold() const {
            return (const T*&)this->naked;
        }
    private:
        ptr    ghost;
        T*     value;
        void*  naked;
    };

}

#endif
