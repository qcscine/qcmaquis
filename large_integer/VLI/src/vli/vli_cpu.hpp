//
//  vli_vpu.hpp
//  vli
//
//  Created by Tim Ewart on 22.03.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

// C - constructors, copy-swap, access operators
template<typename BaseInt, std::size_t Size>
vli_cpu<BaseInt, Size>::vli_cpu(){
    for(size_type i=0; i<size; ++i)
        data_[i] = 0;
}

template<typename BaseInt, std::size_t Size>
vli_cpu<BaseInt, Size>::vli_cpu(int num) {            
 //   assert( static_cast<BaseInt>((num<0) ? -num : num)  < static_cast<BaseInt>(max_value<BaseInt>::value) );
    data_[0] = num; 
    int sign = 0x01 & (num>>(sizeof(int)*8-1));
    for(size_type i=1; i<size; ++i){
        data_[i] = sign*(base<BaseInt>::value+data_mask<BaseInt>::value);
    }
}

template<typename BaseInt, std::size_t Size>
vli_cpu<BaseInt, Size>::vli_cpu(vli_cpu const& r){
    for(size_type i=0; i<size; ++i)
        data_[i] = r.data_[i];
}

template<typename BaseInt, std::size_t Size>
void vli_cpu<BaseInt, Size>::swap(vli_cpu& a, vli_cpu& b){
    boost::swap(a.data_,b.data_); // gcc 4.2, std::swap does not work with static array
}

template<typename BaseInt, std::size_t Size>
vli_cpu<BaseInt, Size>& vli_cpu<BaseInt, Size>::operator= (vli_cpu r){
    swap(*this,r);
    return *this;
}

template<typename BaseInt, std::size_t Size>
BaseInt& vli_cpu<BaseInt, Size>::operator[](size_type i){
    assert( i < size );
    return *(data_+i);
}

template<typename BaseInt, std::size_t Size>
const BaseInt& vli_cpu<BaseInt, Size>::operator[](size_type i) const{
    assert( i < size );
    return *(data_+i);
}
// c equality operator

template<typename BaseInt, std::size_t Size>
vli_cpu<BaseInt, Size> vli_cpu<BaseInt, Size>::operator-() const{
    vli_cpu tmp(*this);
    tmp.negate();
    return tmp;
}


template<typename BaseInt, std::size_t Size>
bool vli_cpu<BaseInt, Size>::operator == (vli_cpu const& vli) const{
    int n = memcmp((void*)data_,(void*)vli.data_,Size*sizeof(BaseInt));
    return (0 == n);
}
        
template<typename BaseInt, std::size_t Size>
bool vli_cpu<BaseInt, Size>::operator != (vli_cpu const& vli) const{
    for(size_type i(0); i < size-1 ; ++i){
        if((*this)[i] != vli[i])
            return true;
    }        
    return false;
}

template<typename BaseInt, std::size_t Size>
bool vli_cpu<BaseInt, Size>::operator < (vli_cpu const& vli) const{
    vli_cpu tmp(*this);
    return ( (tmp-=vli).is_negative() );
}

template<typename BaseInt, std::size_t Size>
bool vli_cpu<BaseInt, Size>::operator < (int i) const{
    vli_cpu tmp1(*this);
    vli_cpu tmp2(i);    
    return ( (tmp1-=tmp2).is_negative() );
}

template<typename BaseInt, std::size_t Size>      
bool vli_cpu<BaseInt, Size>::operator > (vli_cpu vli) const{
    return ( (vli-=*this).is_negative() );
}

template<typename BaseInt, std::size_t Size>
bool vli_cpu<BaseInt, Size>::operator > (int i) const{
    vli_cpu tmp(i);
    return ( (tmp-=*this).is_negative() );
}

// c - negative number

template<typename BaseInt, std::size_t Size>
void vli_cpu<BaseInt, Size>::negate(){
    for(size_type i=0; i < size-1; ++i)
        data_[i] = (~data_[i]);
    data_[size-1] = (~data_[size-1])&(base<BaseInt>::value+data_mask<BaseInt>::value);
    (*this)+=1;
}

template<typename BaseInt, std::size_t Size>
bool vli_cpu<BaseInt, Size>::is_negative() const{
    return static_cast<bool>((data_[size-1]>>data_bits<BaseInt>::value));
}
// c - basic operators
template<typename BaseInt, std::size_t Size>
vli_cpu<BaseInt, Size>& vli_cpu<BaseInt, Size>::operator += (vli_cpu<BaseInt, Size> const& vli){
    using vli::plus_assign;
    plus_assign(*this,vli);
    return *this;
}

template<typename BaseInt, std::size_t Size>
vli_cpu<BaseInt, Size>& vli_cpu<BaseInt, Size>::operator += (BaseInt const a){
    using vli::plus_assign;
    plus_assign(*this,a);
    return *this;
}

template<typename BaseInt, std::size_t Size>
vli_cpu<BaseInt, Size>& vli_cpu<BaseInt, Size>::operator -= (vli_cpu<BaseInt, Size> const& vli){
    using vli::minus_assign;
    minus_assign(*this,vli);
    return *this;
}

template<typename BaseInt, std::size_t Size>
vli_cpu<BaseInt, Size>& vli_cpu<BaseInt, Size>::operator -= (BaseInt const a){
    using vli::minus_assign;
    minus_assign(*this,a);
    return *this;
}

template<typename BaseInt, std::size_t Size>
vli_cpu<BaseInt, Size>& vli_cpu<BaseInt, Size>::operator *= (BaseInt const a){
    using vli::multiplies_assign;
    multiplies_assign(*this,a);
    return *this;
}

template<typename BaseInt, std::size_t Size>
vli_cpu<BaseInt, Size>& vli_cpu<BaseInt, Size>::operator *= (vli_cpu<BaseInt, Size> const& vli){
    using vli::multiplies_assign;
    multiplies_assign(*this,vli);
    return *this;
}

template<typename BaseInt, std::size_t Size>
void vli_cpu<BaseInt, Size>::print_raw(std::ostream& os) const{
    os << "(" ;
    for(size_type i = size-1; i > 0; --i)
        os << std::showbase << std::hex << data_[i]<<" ";
    os << std::showbase << std::hex << data_[0];
    os << ")"<<std::dec;
}    

template<typename BaseInt, std::size_t Size>
void vli_cpu<BaseInt, Size>::print(std::ostream& os) const{
    os<<get_str();
}

/**
 * Returns a string with a base10 represenation of the VLI
 */
template<typename BaseInt, std::size_t Size>
std::string vli_cpu<BaseInt, Size>::get_str() const {
    vli_cpu<BaseInt, size> tmp;
    
    if((*this).is_negative()){
        const_cast<vli_cpu<BaseInt, Size> & >(*this).negate();

        for(size_type i=0; i<size; ++i)
            tmp[i] = (*this)[i];
    
        tmp.negate();
        const_cast<vli_cpu<BaseInt, Size> & >(*this).negate();
    }else{
        for(size_type i=0; i<size; ++i)
            tmp[i] = (*this)[i];
    }
    
    if(tmp.is_negative()){
        tmp.negate();
        size_type ten_exp = order_of_magnitude_base10(tmp);
        return std::string("-")+get_str_helper_inplace(tmp,ten_exp);
    }else{
        size_type ten_exp = order_of_magnitude_base10(tmp);
        return get_str_helper_inplace(tmp,ten_exp);
    }
}

/**
 * Returns the order of magnitude of 'value' in base10
 */
template<typename BaseInt, std::size_t Size>
typename vli_cpu<BaseInt, Size>::size_type vli_cpu<BaseInt, Size>::order_of_magnitude_base10(vli_cpu<BaseInt,size> const& value) const {
    assert(!value.is_negative());
    
    vli_cpu<BaseInt,size> value_cpy(value);
    vli_cpu<BaseInt,size> decimal(1);
    size_type exp = 0;
    
    // Find correct order (10^exp) 
    while(!value_cpy.is_negative()){
        value_cpy=value; // reset
        vli_cpu<BaseInt,size> previous_decimal(decimal);
        decimal *= 10;
        ++exp;
        if(decimal < previous_decimal) // Overflow! (we can handle it.)
        {
            break;
        }
        value_cpy-=decimal;
    }
    --exp;
    return exp;
}

/**
 * A helper function to generate the base10 representation for get_str().
 */
template<typename BaseInt, std::size_t Size>
std::string vli_cpu<BaseInt, Size>::get_str_helper_inplace(vli_cpu<BaseInt,size>& value, size_type ten_exp) const {
    assert(!value.is_negative());
    
    // Create a number 10^(exponent-1) sin
    vli_cpu<BaseInt,size> dec(1);
    for(size_type e=0; e < ten_exp; ++e)
        dec *= 10;
    
    // Find the right digit for 10^ten_exp
    vli_cpu<BaseInt,size> value_cpy(value);
    int digit=0;
    while((!value_cpy.is_negative()) && digit<=11){
        value_cpy = value; // reset
        ++digit;
        if(digit*dec < (digit-1)*dec){ // Overflow (we can handle it.)
            break;
        }
        value_cpy-= digit*dec;
    }
    --digit; // we went to far
    
    assert(digit >=0);
    assert(digit < 10); 
    
    value-= digit*dec;
    
    if(ten_exp <= 0)
        return boost::lexical_cast<std::string>(digit);
    else
        return boost::lexical_cast<std::string>(digit)+get_str_helper_inplace(value,ten_exp-1);
}

// free function algebra 

template <class BaseInt, std::size_t Size>
void mul(vli_cpu<BaseInt, 2*Size>& vli_res, vli_cpu<BaseInt, Size> const&  vli_a, vli_cpu<BaseInt, Size> const& vli_b) { 
    multiplies<BaseInt, Size>(vli_res, vli_a, vli_b);
}

template <class BaseInt, std::size_t Size>
void muladd(vli_cpu<BaseInt, 2*Size>& vli_res, vli_cpu<BaseInt, Size> const&  vli_a, vli_cpu<BaseInt, Size> const& vli_b) {         
        multiply_add_assign<BaseInt, Size>(vli_res, vli_a, vli_b);
}

template <class BaseInt, std::size_t Size>
const vli_cpu<BaseInt, Size> operator + (vli_cpu<BaseInt, Size> vli_a, vli_cpu<BaseInt, Size> const& vli_b){
    vli_a += vli_b;
    return vli_a;
}
    
template <class BaseInt, std::size_t Size>
const vli_cpu<BaseInt, Size> operator + (vli_cpu<BaseInt, Size> vli_a, int b){
    vli_a += b;
    return vli_a;
}
    
template <class BaseInt, std::size_t Size>
const vli_cpu<BaseInt, Size> operator + (int b, vli_cpu<BaseInt, Size> const& vli_a){
    return vli_a+b;
}
    
template <class BaseInt, std::size_t Size>
const vli_cpu<BaseInt, Size> operator - (vli_cpu<BaseInt, Size> vli_a, vli_cpu<BaseInt, Size> const& vli_b){
    vli_a -= vli_b;
    return vli_a;
}
    
template <class BaseInt, std::size_t Size>
const vli_cpu<BaseInt, Size> operator - (vli_cpu<BaseInt, Size> vli_a, int b){
    vli_a -= b;
    return vli_a;
}

template <class BaseInt, std::size_t Size>
const vli_cpu<BaseInt, Size> operator * (vli_cpu<BaseInt, Size>  vli_a, vli_cpu<BaseInt, Size> const& vli_b){
    vli_a *= vli_b;
    return vli_a;
}

template <class BaseInt, std::size_t Size>
const vli_cpu<BaseInt, Size> operator * (vli_cpu<BaseInt, Size> vli_a, int b){
    vli_a *= b;
    return vli_a;
}

template <class BaseInt, std::size_t Size>
const vli_cpu<BaseInt, Size> operator * (int b, vli_cpu<BaseInt, Size> const& a){
    return a*b;
}

//stream
template<typename BaseInt, std::size_t Size>
std::ostream& operator<< (std::ostream& os,  vli_cpu<BaseInt, Size> const& vli){
    vli.print_raw(os);
    return os;
}





