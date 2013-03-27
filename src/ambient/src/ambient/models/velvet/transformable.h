#ifndef AMBIENT_MODELS_VELVET_TRANSFORMABLE
#define AMBIENT_MODELS_VELVET_TRANSFORMABLE

namespace ambient { namespace models { namespace velvet {

    class transformable {
    public:
        union numeric_union { 
            bool b; 
            double d; 
            std::complex<double> c; 
            operator bool& ();
            operator double& ();
            operator std::complex<double>& ();
            void operator = (bool value);
            void operator = (double value);
            void operator = (std::complex<double> value);
            numeric_union(){ }
        };
        void* operator new  (size_t, void*);
        void operator delete (void*, void*);
        virtual numeric_union eval() const = 0;
        virtual transformable& operator += (transformable& r) = 0;

        const transformable* l;
        const transformable* r;
        mutable numeric_union v;
        void* generator;
        int sid;
    };

    template <typename T>
    class transformable_value : public transformable {
    public:
        transformable_value(T value);
        virtual transformable::numeric_union eval() const;
        virtual transformable& operator += (transformable& r);
    };

    template <typename T, typename FP, FP OP>
    class transformable_expr : public transformable {
    public:
        transformable_expr(const transformable* l);
        transformable_expr(const transformable* l, const transformable* r);
        virtual transformable::numeric_union eval() const;
        virtual transformable& operator += (transformable& r); 
    };

} } }

#endif
