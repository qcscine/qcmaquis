#ifndef VLI_MAIN_MISC_HPP
#define VLI_MAIN_MISC_HPP
namespace vli
{
template <typename Polynomial>
struct polynomial_multiply_type_gmp {
};

template <typename Coeff, int Order, class Var0, class Var1, class Var2, class Var3>
struct polynomial_multiply_type_gmp<polynomial<Coeff,max_order_each<Order>,Var0,Var1,Var2,Var3> > {
    typedef polynomial<mpz_class,max_order_each<Order>,Var0,Var1,Var2,Var3> type;
    typedef polynomial<mpz_class,max_order_each<2*Order>,Var0,Var1,Var2,Var3> type_res;
};

template <typename Coeff, int Order, class Var0, class Var1, class Var2, class Var3>
struct polynomial_multiply_type_gmp<polynomial<Coeff,max_order_combined<Order>,Var0,Var1,Var2,Var3> > {
    typedef polynomial<mpz_class,max_order_combined<Order>,Var0,Var1,Var2,Var3> type;
    typedef polynomial<mpz_class,max_order_combined<2*Order>,Var0,Var1,Var2,Var3> type_res;
};
}

namespace tools
{
    template <typename Polynomial>
    void converter(typename vli::vector<Polynomial>& v_vli,  vli::vector<typename vli::polynomial_multiply_type_gmp<Polynomial>::type>& v_gmp ){
        typedef typename vli::polynomial_multiply_type_gmp<Polynomial>::type Polynomial_gmp;
        typename Polynomial::iterator it_poly_vli;
        typename Polynomial_gmp::iterator it_poly_gmp;
    
        std::size_t vec_size = v_vli.size();
    
        for(int i=0; i< vec_size; ++i){
            it_poly_vli = v_vli[i].begin();
            it_poly_gmp = v_gmp[i].begin();
            for(; it_poly_vli != v_vli[i].end(); ++it_poly_vli, ++it_poly_gmp){
                (*it_poly_gmp) = vli::detail::gmp_convert_helper<mpz_class>::apply(*it_poly_vli);
            }
        }
    }

    template <typename Polynomial>
    bool equal(typename vli::polynomial_multiply_result_type<Polynomial>::type& p_vli, typename vli::polynomial_multiply_type_gmp<Polynomial>::type_res& p_gmp){
        typedef typename vli::polynomial_multiply_type_gmp<Polynomial>::type Polynomial_gmp;
        typedef typename vli::polynomial_multiply_type_gmp<Polynomial>::type_res Polynomial_gmp_res;

        typename vli::polynomial_multiply_result_type<Polynomial>::type::iterator it_poly_vli = p_vli.begin(); 

        typename vli::polynomial_multiply_type_gmp<Polynomial>::type_res p_gmp_tmp;
        typename vli::polynomial_multiply_type_gmp<Polynomial>::type_res::iterator it_poly_gmp = p_gmp_tmp.begin(); 

        for(; it_poly_vli != p_vli.end(); ++it_poly_vli, ++it_poly_gmp)
            (*it_poly_gmp) = vli::detail::gmp_convert_helper<mpz_class>::apply(*it_poly_vli);

        return (p_gmp == p_gmp_tmp);
    }



}
#endif // VLI_MAIN_MISC_HPP
