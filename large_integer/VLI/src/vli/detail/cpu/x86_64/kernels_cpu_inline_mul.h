
//  inline_mul.h
//  vli
//
//  Created by Timothée Ewart on 22.01.13.
//
//

#ifndef vli_inline_mul_h
#define vli_inline_mul_h

namespace vli{
    namespace detail{

// WARNING : ONLY work with signed integer

/* \cond I do not need this part in the doc*/

template<int num>
inline void inline_helper_mul(boost::uint64_t* x,  const boost::uint64_t  a, const boost::uint64_t  b);

template<>
inline void inline_helper_mul<5>(boost::uint64_t* x,  const boost::uint64_t  a, const boost::uint64_t  b){
    boost::uint64_t tmp_rax, tmp_rdx; //fictive registers
    asm __volatile__(
                     "mulq %4 \n\t" // b %%rdx
                     "addq %[tmp_rax], (%[x])   \n\t"
                     "adcq %[tmp_rdx], 8(%[x])  \n\t"
                     "adcq $0        , 16(%[x]) \n\t"
                     "adcq $0        , 24(%[x]) \n\t"
                     "adcq $0        , 32(%[x]) \n\t"
                     : [tmp_rax] "=&a"(tmp_rax), [tmp_rdx] "=&d"(tmp_rdx)
                     : [x] "r" (x), "%0" (a), "r" (b) /* a into rax directly */
                     : "memory", "cc");
}

template<>
inline void inline_helper_mul<4>(boost::uint64_t* x,  const boost::uint64_t  a, const boost::uint64_t  b){
    boost::uint64_t tmp_rax, tmp_rdx; //fictive registers
    asm __volatile__(
                     "mulq %4 \n\t" // b %%rdx
                     "addq %[tmp_rax], (%[x])   \n\t"
                     "adcq %[tmp_rdx], 8(%[x])  \n\t"
                     "adcq $0        , 16(%[x]) \n\t"
                     "adcq $0        , 24(%[x]) \n\t"
                     : [tmp_rax] "=&a"(tmp_rax), [tmp_rdx] "=&d"(tmp_rdx)
                     : [x] "r" (x), "%0" (a), "r" (b) /* a into rax directly */
                     : "memory", "cc");
}

template<>
inline void inline_helper_mul<3>(boost::uint64_t* x,  const boost::uint64_t  a, const boost::uint64_t  b){
    boost::uint64_t tmp_rax, tmp_rdx; //fictive registers
    asm __volatile__(
                     "mulq %4 \n\t"
                     "addq %[tmp_rax], (%[x])   \n\t"
                     "adcq %[tmp_rdx], 8(%[x])  \n\t"
                     "adcq $0        , 16(%[x]) \n\t"
                     : [tmp_rax] "=&a"(tmp_rax), [tmp_rdx] "=&d"(tmp_rdx)
                     : [x] "r" (x), "%0" (a), "r" (b) /* a into rax directly */
                     : "memory", "cc");
}

template<>
inline void inline_helper_mul<2>(boost::uint64_t* x,  const boost::uint64_t  a, const boost::uint64_t  b){
    boost::uint64_t tmp_rax, tmp_rdx; //fictive registers
    asm __volatile__(
                     "mulq %4 \n\t"
                     "addq  %[tmp_rax], (%[x]) \n\t"
                     "adcq  %[tmp_rdx], 8(%[x]) \n\t"
                     : [tmp_rax] "=&a"(tmp_rax), [tmp_rdx] "=&d"(tmp_rdx)
                     : [x] "r" (x), [a] "%0" (a), [b] "r" (b) /* a into rax directly */
                     : "memory", "cc");
}

template<std::size_t numbits> // 128 -> 256
inline void inline_mul_extend(boost::uint64_t* z,  const boost::uint64_t *x, const boost::uint64_t *y);


template<> // 128 -> 256
inline void inline_mul_extend<2>(boost::uint64_t* z,  const boost::uint64_t *x, const boost::uint64_t *y){
    //clean z to avoid pb
    z[0] ^= z[0];
    z[1] ^= z[1];
    z[2] ^= z[2];
    z[3] ^= z[3];
    
    inline_helper_mul<2>(&z[0],x[0],y[0]);
    inline_helper_mul<2>(&z[1],x[0],y[1]);
    
    inline_helper_mul<3>(&z[1],x[1],y[0]);
    inline_helper_mul<2>(&z[2],x[1],y[1]);
}

template<> // 192 -> 384
inline void inline_mul_extend<3>(boost::uint64_t* z,  const boost::uint64_t *x, const boost::uint64_t *y){
    //clean z to avoid pb
    z[0] ^= z[0];
    z[1] ^= z[1];
    z[2] ^= z[2];
    z[3] ^= z[3];
    z[4] ^= z[4];
    z[5] ^= z[5];
    
    inline_helper_mul<2>(&z[0],x[0],y[0]);
    inline_helper_mul<2>(&z[1],x[0],y[1]);
    inline_helper_mul<2>(&z[2],x[0],y[2]);
    
    inline_helper_mul<2>(&z[3],x[1],y[2]);
    inline_helper_mul<3>(&z[1],x[1],y[0]);
    inline_helper_mul<3>(&z[2],x[1],y[1]);
    
    inline_helper_mul<2>(&z[4],x[2],y[2]);
    inline_helper_mul<3>(&z[2],x[2],y[0]);
    inline_helper_mul<3>(&z[3],x[2],y[1]);
}


template<> // 256 -> 512
inline void inline_mul_extend<4>(boost::uint64_t* z,  const boost::uint64_t *x, const boost::uint64_t *y){
    //clean z to avoid pb
    z[0] ^= z[0];
    z[1] ^= z[1];
    z[2] ^= z[2];
    z[3] ^= z[3];
    z[4] ^= z[4];
    z[5] ^= z[5];
    z[6] ^= z[6];
    z[7] ^= z[7];
    
    
    inline_helper_mul<2>(&z[0],x[0],y[0]);
    inline_helper_mul<2>(&z[1],x[0],y[1]);
    inline_helper_mul<2>(&z[2],x[0],y[2]);
    inline_helper_mul<2>(&z[3],x[0],y[3]);
    
    inline_helper_mul<2>(&z[4],x[1],y[3]);
    inline_helper_mul<3>(&z[3],x[1],y[2]);
    inline_helper_mul<3>(&z[2],x[1],y[1]);
    inline_helper_mul<4>(&z[1],x[1],y[0]);
    
    inline_helper_mul<2>(&z[5],x[2],y[3]);
    inline_helper_mul<3>(&z[4],x[2],y[2]);
    inline_helper_mul<3>(&z[3],x[2],y[1]);
    inline_helper_mul<4>(&z[2],x[2],y[0]);
    
    inline_helper_mul<2>(&z[6],x[3],y[3]);
    inline_helper_mul<3>(&z[5],x[3],y[2]);
    inline_helper_mul<3>(&z[4],x[3],y[1]);
    inline_helper_mul<4>(&z[3],x[3],y[0]);
}

/* \endcond */

    }} //end namespace
 
#endif
