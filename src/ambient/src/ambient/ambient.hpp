#ifndef AMBIENT
#define AMBIENT
// {{{ system includes
#include <mpi.h>
#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string>
#include <limits>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <memory.h>
#include <stdarg.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <algorithm>
#include <pthread.h>
// #include <execinfo.h>
// }}}

#ifdef AMBIENT_CILK
    #include <cilk/cilk.h>
    #define AMBIENT_NUM_THREADS __cilkrts_get_total_workers()
    #define AMBIENT_THREAD_ID __cilkrts_get_worker_number()
    #define AMBIENT_THREAD cilk_spawn
    #define AMBIENT_SMP_ENABLE
    #define AMBIENT_SMP_DISABLE
#elif defined(AMBIENT_OMP)
    #include <omp.h>
    #define AMBIENT_THREAD_ID omp_get_thread_num()
    #define AMBIENT_PRAGMA(a) _Pragma( #a )
    #define AMBIENT_THREAD AMBIENT_PRAGMA(omp task untied)
    #define AMBIENT_SMP_ENABLE AMBIENT_PRAGMA(omp parallel) { AMBIENT_PRAGMA(omp single)
    #define AMBIENT_SMP_DISABLE }
    #define AMBIENT_NUM_THREADS [&]()->int{ int n; AMBIENT_SMP_ENABLE \
                                { n = omp_get_num_threads(); } \
                                AMBIENT_SMP_DISABLE return n; }()
#else
    #define AMBIENT_NUM_THREADS 1
    #define AMBIENT_THREAD_ID   0
    #define AMBIENT_THREAD
    #define AMBIENT_SMP_ENABLE
    #define AMBIENT_SMP_DISABLE
#endif

//#define AMBIENT_COMPUTATIONAL_TIMINGS
//#define AMBIENT_CHECK_BOUNDARIES
//#define AMBIENT_TRACE void* b[10]; backtrace_symbols_fd(b,backtrace(b,10),2);
//#define AMBIENT_LOOSE_FUTURE

#define AMBIENT_NUM_PROCS             2
#define AMBIENT_MAX_SID               2147483647
#define AMBIENT_STACK_RESERVE         65536
#define AMBIENT_COLLECTOR_STR_RESERVE 65536
#define AMBIENT_COLLECTOR_RAW_RESERVE 1024
#define AMBIENT_SCOPE_SWITCH_FACTOR   20480
#define AMBIENT_BULK_CHUNK            41943040
#define AMBIENT_FUTURE_SIZE           64
#define AMBIENT_IB                    512

#define BULK_REGION      0
#define DELEGATED_REGION 0 // same as bulk - don't deallocate
#define DEFAULT_REGION   1
#define PERSIST_REGION   13

namespace ambient {
    inline int get_num_threads(){
        static int n = AMBIENT_NUM_THREADS; return n;
    }
    enum complexity { N, N2, N3 };
    enum locality   { remote, local, common };
    enum scope_t    { base, single, shared  };
    enum memory_t   { staged, normal, leak  };
}

#include "ambient/utils/memory.hpp"
#include "ambient/models/velvet/model.h"
#include "ambient/channels/mpi/channel.h"
#include "ambient/controllers/velvet/controller.h"
#include "ambient/utils/auxiliary.hpp"
#include "ambient/utils/io.hpp"
#include "ambient/interface/typed.hpp"
#include "ambient/interface/kernel.hpp"
#include "ambient/interface/access.hpp"
#endif
