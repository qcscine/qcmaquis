#include "ambient/ambient.hpp"

template<typename T> class vector;

namespace detail { using namespace ambient;

    template<typename T>
    void init_value(unbound< vector<T> >& a, T& value){
        size_t size = get_length(a)-1;
        T* a_ = versioned(a).elements;
        for(size_t i = 0; i < size; ++i) a_[i] = value;
        versioned(a).count = 0;
    }

    template<typename T>
    void add(vector<T>& a, const vector<T>& b){
        size_t size = get_length(a)-1;
        T* a_ = versioned(a).elements;
        T* b_ = versioned(b).elements;
        for(size_t i = 0; i < size; ++i) a_[i] += b_[i];
        versioned(a).count = versioned(a).count+1;
    }
}

ambient_reg(detail::init_value, init_value)
ambient_reg(detail::add, add)

template <typename T>
class vector {
public:
    class iterator {
    public:
        typedef T value_type;
        iterator(vector<T>& base) : base(base) {}
        vector<T>& base;
    };

    ambient_version(
        T count;
        T elements[ AMBIENT_VAR_LENGTH ];
    );
    vector(size_t length, T value) : versioned(length+1, sizeof(T)) {
        init_value<T>::spawn(*this, value);
    }
    iterator begin(){
        return iterator(*this);
    }
    iterator end(){
        return iterator(*this);
    }
};

namespace ambient {

    template<class InputIterator, class Function>
    void for_each (InputIterator first, InputIterator last, Function fn){
        ambient::lambda([fn](vector<int>& a){ 
            int* a_ = versioned(a).elements;
            std::for_each(a_, a_+get_length(a)-1, fn);
        })(first.base);
    }

    template <class InputIterator, class OutputIterator, class BinaryOperation>
    void transform (InputIterator first1, InputIterator last1,
                    InputIterator first2, 
                    OutputIterator result,
                    BinaryOperation binary_op)
    {
        ambient::lambda([binary_op](const vector<int>& first1_, const vector<int>& first2_, unbound< vector<int> >& result_){
            int* ifirst1 = versioned(first1_).elements;
            int* ifirst2 = versioned(first2_).elements;
            int* iresult = versioned(result_).elements;
            std::transform(ifirst1, ifirst1+get_length(first1_)-1, ifirst2, iresult, binary_op);
        })(first1.base, first2.base, result.base);
    }

    template <class InputIterator, class OutputIterator, class UnaryOperation>
    void transform (InputIterator first1, InputIterator last1,
                    OutputIterator result, UnaryOperation op)
    {
        ambient::lambda([op](const vector<int>& first1_, unbound< vector<int> >& result_){
            int* ifirst1 = versioned(first1_).elements;
            int* iresult = versioned(result_).elements;
            std::transform(ifirst1, ifirst1+get_length(first1_)-1, iresult, op);
        })(first1.base, result.base);
    }

    template<class InputIterator>
    void sequence (InputIterator first, InputIterator last){
    }

    template <class ForwardIterator, class T>
    void fill (ForwardIterator first, ForwardIterator last, const T& val){
    }

    template <class ForwardIterator, class T>
    void replace (ForwardIterator first, ForwardIterator last, const T& old_value, const T& new_value){
    }

    template <class ForwardIterator, class Generator>
    void generate (ForwardIterator first, ForwardIterator last, Generator gen){
    }

    template <class InputIterator, class OutputIterator>
    OutputIterator copy (InputIterator first, InputIterator last, OutputIterator result){
    }

    template <class ForwardIterator, class T>
    ForwardIterator remove (ForwardIterator first, ForwardIterator last, const T& val){
    }

    template <class ForwardIterator>
    ForwardIterator unique (ForwardIterator first, ForwardIterator last){
    }

    template <class ForwardIterator, class BinaryPredicate>
    ForwardIterator unique (ForwardIterator first, ForwardIterator last, BinaryPredicate pred){
    }

    template <class RandomAccessIterator>
    void sort (RandomAccessIterator first, RandomAccessIterator last){
    }

    template <class RandomAccessIterator, class Compare>
    void sort (RandomAccessIterator first, RandomAccessIterator last, Compare comp){
    }

    template< class InputIt, class T>
    T reduce( InputIt first, InputIt last, T init ){
    }

    template< class InputIt, class T, class BinaryOperation>
    T reduce( InputIt first, InputIt last, T init, BinaryOperation op ){
    }


    template <class InputIterator, class T>
    InputIterator find (InputIterator first, InputIterator last, const T& val){
    }

}


int main(){ using namespace ambient;

    vector<int> a(10, 13);
    vector<int> b(10, 10);
                                 { scope<> select(1);
    add<int>::spawn(a, b);
                                 }
                                 { scope<> select(0); 
    add<int>::spawn(a, b);
                                 }
    for(int i = 0; i < 10; i++)
    cout << "After sync: " << get(a).elements[i] << "; count: " << get(a).count << "\n";

    cout << "Messing with lambda:\n";
    int num = 13;

    ambient::for_each( a.begin(), a.end(), [&] (int& val){ val += num; } );
    ambient::for_each( a.begin(), a.end(), [&] (int& val){ val += num; } );

    ambient::transform( a.begin(), a.end(), b.begin(), a.begin(), [](double a, double b){ return a + b; } );

    ambient::lambda([&](const vector<int>& val, const vector<int>& val2){ 
            size_t size = get_length(val)-1;
            int* val_ = versioned(val).elements;
            int* val2_ = versioned(val2).elements;
            for(size_t i = 0; i < size; ++i) std::cout << val_[i] << " vs " << val2_[i] << "; num is " << num << "\n";
    })(a, b);

    ambient::lambda([&](const vector<int>& val){ 
            size_t size = get_length(val)-1;
            int* val_ = versioned(val).elements;
            for(size_t i = 0; i < size; ++i) std::cout << val_[i] << "\n";
    })(a);


    ambient::sync();

    /*// possible todo:
    ambient::for_(a.begin(), b.begin(), [&](parallel_iterator& a_, parallel_iterator& b_)
    { 
        for(int i = 0; i < get_length(a); i++){
            a_[i] += b_[i];
        }
    });*/

    return 0;
}

