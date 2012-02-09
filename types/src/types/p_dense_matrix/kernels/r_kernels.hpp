namespace ambient {

    template <typename T>
    void plus_reduce(void* dst, void* src){
        assert(false); // only partially specialized
    }

    template<>
    void plus_reduce< maquis::types::p_dense_matrix_impl<double> >(void* dst, void* src){
        double* a = (double*)dst;
        double* u = (double*)src;
        *a += *u; // for now - atomic
    }

}
