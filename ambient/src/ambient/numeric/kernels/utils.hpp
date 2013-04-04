#ifndef AMBIENT_NUMERIC_MATRIX_KERNELS_UTILS
#define AMBIENT_NUMERIC_MATRIX_KERNELS_UTILS

namespace ambient {

    template<typename T> inline dim2 dim(T& ref){ 
        return ref.core->spec.dim;
    }

    template<typename T> inline size_t square_dim(T& ref){ 
        return ref.core->spec.dim.square();
    }

    template<typename T> inline size_t size(T& ref){ 
        return ref.core->spec.size;  
    }

    template <typename T>
    inline void memcpy(T* dd, T *sd, size_t w, T alfa){
        std::memcpy(dd, sd, w);
    }

    template <typename T>
    inline void memscal(T* dd, T *sd, size_t w, T alfa){
        int z = w/sizeof(T);
        do{ *dd++ = alfa*(*sd++); }while(--z > 0); // note: dd != sd
    }

    template <typename T>
    inline void memscala(T* dd, T *sd, size_t w, T alfa){
        int z = w/sizeof(T);
        do{ *dd++ += alfa*(*sd++); }while(--z > 0); // note: dd != sd
    }

    template <typename T> inline T dot(T* a, T* b, int size){
        T summ(0);
        for(size_t k=0; k < size; k++)
           summ += a[k]*b[k];
        return summ;
    }

    inline double dot(double* a, double* b, int size){
        static const int ONE = 1;
        return ddot_(&size, a, &ONE, b, &ONE);
    }

    template<typename T, void(*PTF)(T* dd, T* sd, size_t w, T alfa)>
    inline void memptf(T* dst, int ldb, dim2 dst_p, 
                       T* src, int lda, dim2 src_p, 
                       dim2 size, T alfa = 0.0)
    {
        #ifdef AMBIENT_CHECK_BOUNDARIES
        if(ambient::dim(dst).x - dst_p.x < size.x || ambient::dim(dst).y - dst_p.y < size.y ||
           ambient::dim(src).x - src_p.x < size.x || ambient::dim(src).y - src_p.y < size.y){
            ambient::cout << "Error: invalid memory movement: \n";
            ambient::cout << "Matrix dst " << ambient::dim(dst).x << "x" << ambient::dim(dst).y << "\n";
            ambient::cout << "Dest p " << dst_p.x << "x" << dst_p.y << "\n";
            ambient::cout << "Matrix src " << ambient::dim(src).x << "x" << ambient::dim(src).y << "\n";
            ambient::cout << "Src p " << src_p.x << "x" << src_p.y << "\n";
            ambient::cout << "Block size " << size.x << "x" << size.y << "\n";
            AMBIENT_TRACE
        }
        #endif
        int n = size.x;
        int m = size.y*sizeof(T);
        T* sd = src + src_p.y + src_p.x*lda;
        T* dd = dst + dst_p.y + dst_p.x*ldb;
        do{ PTF(dd, sd, m, alfa); sd += lda; dd += ldb; }while(--n > 0);
    }
}

#endif
