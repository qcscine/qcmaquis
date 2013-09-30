#ifndef PARALLEL_FOR_HPP
#define PARALLEL_FOR_HPP

#ifdef AMBIENT
    typedef ambient::scope<ambient::single> locale;
    #define parallel_for(constraint, ...) constraint; for(__VA_ARGS__)
    #define semi_parallel_for(constraint, ...) constraint; for(__VA_ARGS__)
#elif defined(MAQUIS_OPENMP)
    typedef std::size_t locale;
    #define parallel_pragma(a) _Pragma( #a )
    #define parallel_for(constraint, ...) parallel_pragma(omp parallel for schedule(dynamic, 1)) for(__VA_ARGS__)
    #define semi_parallel_for(constraint, ...) for(__VA_ARGS__)
#else
    typedef std::size_t locale;
    #define parallel_for(constraint, ...) for(__VA_ARGS__)
    #define semi_parallel_for(constraint, ...) for(__VA_ARGS__)
#endif

#endif
