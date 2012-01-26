#ifndef VLI_ALGORITHMS_POLYNOMIAL_CPU
#define VLI_ALGORITHMS_POLYNOMIAL_CPU

namespace vli{ 

template<class Vli, unsigned int Order>
class polynomial_cpu;
    
template<class BaseInt, std::size_t Size>
class vli_cpu;
    
    
/** Very important ! All these algo used the truncation version **/
  
/** First Algo based on block : n threads possible **/    
    
template <class BaseInt, std::size_t Size, unsigned int Order>
void triangle_up(unsigned int block_ai, int block_bj,
                 polynomial_cpu<vli_cpu<BaseInt, Size>, 2*Order> & result,    
                 polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p1, 
                 polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p2)
{
    unsigned int n(0);
    unsigned int offset_block_ai = block_ai*Order;
    unsigned int offset_block_bj = block_bj*Order;
    
    for(unsigned int i = 0; i < Order-1; ++i){
        for(unsigned int j = n; j < Order-1; ++j){
            result.coeffs_[(offset_block_ai+offset_block_bj)*2+j] += p1.coeffs_[j-n+offset_block_ai] * p2.coeffs_[n+offset_block_bj];
        }    
        n++;
    }        
};

template <class BaseInt, std::size_t Size, unsigned int Order>
void diag(unsigned block_ai, unsigned block_bj,
          polynomial_cpu<vli_cpu<BaseInt, Size>, 2*Order> & result,    
          polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p1, 
          polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p2)
{
    unsigned int OrderMinusOne = Order-1;
    unsigned int offset_block_ai = block_ai*Order;
    unsigned int offset_block_bj = block_bj*Order;    
    
    for(unsigned int i = 0; i < Order; ++i)
        result.coeffs_[(offset_block_ai+offset_block_bj)*2+OrderMinusOne] += p1.coeffs_[offset_block_bj+Order-1-i] * p2.coeffs_[offset_block_ai+i];
}

template <class BaseInt, std::size_t Size, unsigned int Order>
void triangle_down(unsigned int block_ai, unsigned int block_bj, 
                   polynomial_cpu<vli_cpu<BaseInt, Size>, 2*Order> & result,    
                   polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p1, 
                   polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p2)
{
    unsigned int n(0);
    unsigned int offset_block_ai = (block_ai+1)*(Order)-1;
    unsigned int offset_block_bj = (block_bj+1)*(Order)-1;
    
    for(unsigned int i = 0; i < Order-1; ++i){
        for(unsigned int j = n; j < Order-1; ++j){
            result.coeffs_[(offset_block_ai+offset_block_bj)*2-2*Order+2-j] +=  p1.coeffs_[offset_block_ai+n-j] * p2.coeffs_[offset_block_bj-n]; // will be VLI *
        }    
        n++;
    }        
}

template <class BaseInt, std::size_t Size, unsigned int Order>
void block_algo(int i, int j, 
                polynomial_cpu<vli_cpu<BaseInt, Size>, 2*Order> & result,    
                polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p1, 
                polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p2)
{
    triangle_up(i,j,result,p1,p2);  
    diag(i,j,result,p1,p2);
    triangle_down(i,j,result,p1,p2);
}

/** second Algo  based on diagonals n*n symetries **/    

template <class BaseInt, std::size_t Size, unsigned int Order>
void diagonal_up(unsigned int n,
                polynomial_cpu<vli_cpu<BaseInt, Size>, 2*Order> & result,    
                polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p1, 
                polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p2)
{
    unsigned int qa,ra,qb,rb,pos; // find all indexes
    
    for(int i(0); i <= n; i++){
        qa = i/Order;
        ra = i%Order;
        qb = (n-i)/Order;
        rb = (n-i)%Order;
        pos = 2*(qa+qb)*Order + (ra+rb);
        result.coeffs_[pos] += p1.coeffs_[n-i]*p2.coeffs_[i];  
    }
}
    
template <class BaseInt, std::size_t Size, unsigned int Order>
void diagonal_down(unsigned int n,
                   polynomial_cpu<vli_cpu<BaseInt, Size>, 2*Order> & result,    
                   polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p1, 
                   polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p2)
{    
    int qa,ra,qb,rb,pos; // find all indexes

    int j = Order*Order-1;
 
    for(int i(Order*Order-n+1); i < Order*Order; i++){
        qa = i/Order;
        ra = i%Order;
        qb = j/Order;
        rb = j%Order;
        pos = 2*(qa+qb)*Order + (ra+rb);
        result.coeffs_[pos] += p1.coeffs_[j]*p2.coeffs_[i];  
        j--;        
    }    
}   
    
}

#endif