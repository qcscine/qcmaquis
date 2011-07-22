
namespace vli
{

template <typename Vli>
struct max_int_value
{
    enum { value = data_mask<typename Vli::value_type>::value };
};

} // namespace vli
