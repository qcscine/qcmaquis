#include <boost/concept_check.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <stdexcept>

namespace blas
{

template <typename X>
struct ResizableMatrix
        : Matrix<X>
{
    public:
    BOOST_CONCEPT_USAGE(ResizableMatrix)
    {
        typename boost::remove_const<X>::type x(1,1);

        // Resize
        resize(x,2,2);

        // Append
        std::vector<typename X::value_type> dataA(2,2);
        std::vector<typename X::value_type> dataB(3,2);
        append_row(x, std::make_pair(dataA.begin(),dataA.end()) );
        append_column(x, std::make_pair(dataB.begin(),dataB.end()) );

        // Remove
        remove_row(x,1);
        remove_column(x,1);

        // Insert
        insert_row(x,1, std::make_pair(dataA.begin(),dataA.end()) );
        insert_column(x,1, std::make_pair(dataB.begin(),dataB.end()) );
    }
};

}
