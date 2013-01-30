#ifndef AMBIENT_UTILS_REDUCE
#define AMBIENT_UTILS_REDUCE

#define __a_ceil(x) (((double)x-(int)x) == 0 ? (int)x : (int)x+1)

inline size_t __a_mod_modern(size_t size, size_t tile){
    size_t mask[2] = {(size & (tile-1)), tile}; 
    return mask[!mask[0]];
}

inline size_t __a_mod(size_t size, size_t tile){
    size_t m = size % tile;
    if(m == 0) m = tile;
    return m;
}

template<typename T>
inline void __a_reduce(std::vector<T*>& seq){
    if(seq.size() == 1) return;
    for(int stride = 1; stride < seq.size(); stride *= 2)
        for(int k = stride; k < seq.size(); k += stride*2){
            *seq[k-stride] += *seq[k];
        }
}

#endif
