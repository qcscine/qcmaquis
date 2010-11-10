#include <general_matrix/general_matrix.hpp>
#include <general_matrix/matrix_interface.hpp>
#include <general_matrix/matrix_iterators.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <CopyOnWrite/aligned_allocator.h>
#include <CopyOnWrite/cow_vector.h>

template <typename T>
class PsudoRandomGenerator
{
    private:
        static boost::mt19937 rng_;
    public:
        T operator()()
        {
            return rng_() %10;
        }
};

template <typename T>
boost::mt19937 PsudoRandomGenerator<T>::rng_(342324);

template <typename T>
bool matrix_matrix_multiplication_test(blas::general_matrix<T> const& a, blas::general_matrix<T> const& b)
{
    blas::general_matrix<T> c( a*b);
    std::cout<<c<<std::endl;
    return false;
}

template <typename T>
bool plus_test(blas::general_matrix<T> const& a, blas::general_matrix<T> const& b)
{
    blas::general_matrix<T> c( a+b);
    std::cout<<c<<std::endl;
    blas::general_matrix<T> d( c-a);
    if( b == d) return true;
    else return false;
}

template <typename T>
void test_iterator_conversions(blas::general_matrix<T> const& a)
{
    std::pair<typename blas::general_matrix<T>::const_column_element_iterator,typename blas::general_matrix<T>::const_column_element_iterator> ita = a.column(0);
    std::pair<typename blas::general_matrix<T>::const_column_element_iterator,typename blas::general_matrix<T>::const_column_element_iterator> itb = a.column(1);
    itb = ita;
    int result =0;
    for(std::size_t i=0; i<a.num_columns(); ++i)
    {
        result = std::accumulate(ita.first, ita.second, result);
    }
    std::cout<<result<<std::endl;

    std::cout<<std::distance(ita.first, ita.second)<<std::endl;
}


int main()
{
    using blas::general_matrix;
    general_matrix<double> a(50,80);
    general_matrix<double> b(80,10);
    general_matrix<double,copy_on_write_vector<double> > c(20,10,0.);
    
    for(std::size_t i=0; i<num_columns(a); ++i)
    {
        std::pair<general_matrix<double>::column_element_iterator,general_matrix<double>::column_element_iterator> range (a.column(i));
        std::generate( range.first, range.second, PsudoRandomGenerator<double>());
    }
    for(std::size_t i=0; i<num_columns(b); ++i)
    {
        std::pair<general_matrix<double>::column_element_iterator,general_matrix<double>::column_element_iterator> range (b.column(i));
        std::generate( range.first, range.second, PsudoRandomGenerator<double>());
    }
    std::cout<<a<<std::endl;
    std::cout<<b<<std::endl;
    matrix_matrix_multiplication_test(a,b);
    a.resize(20,10);
    b.resize(20,10);
    general_matrix<double> d(c);
    plus_test(a,d);
    test_iterator_conversions(a);
    std::cout<<std::boolalpha<<plus_test(a,b)<<std::endl;

    blas::vector<float> e(50,0.1f);
    blas::vector<float,copy_on_write_vector<float> > f(30,0.2f);
    e=f;
    std::cout<<e<<std::endl;
    return 0;
}

