
namespace vli
{

namespace utils
{
template <typename Vector>
class entrywise
{
    public:
        Vector const& vector;

        explicit entrywise(Vector const& v)
            :vector(v)
        {
        }
};
}

template <typename Vector>
inline utils::entrywise<Vector> entrywise(Vector const& v)
{
    return utils::entrywise<Vector>(v);
}

}
