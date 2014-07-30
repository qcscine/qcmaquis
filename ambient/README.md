Ambient environment variables:


AMBIENT_VERBOSE                        # print-out Ambient configuration prior to running
[not set]                              #

AMBIENT_MKL_NUM_THREADS=[thread count] # enable selective threading: MKL will use only 1 thread 
[not set]                              # unless sync is called passing mkl_parallel() - then the
                                       # [thread count] is used.

MKL_NUM_THREADS=[thread count]         # MKL will use [thread count] in all of its calls.
[auto]                                 # Warning: can cause performance degradation due to
                                       # resources overloading. Use AMBIENT_MKL_NUM_THREADS 
                                       # to override it.

AMBIENT_BULK_LIMIT=[chunks number]     # limit the data bulk memory consumption by [chunks number]
[unlimited]                            #

AMBIENT_BULK_REUSE                     # setting this variable will enable bulk garbage collection
[not set]                              #

AMBIENT_BULK_DEALLOCATE                # deallocate data bulk every time the sync has finished
[not set]                              #


Resource overloading impact on a tiny test-case:

115 s - MPI_NUM_PROCS=12                           run dmrg short 
115 s - MPI_NUM_PROCS=12 AMBIENT_MKL_NUM_THREADS=1 run dmrg short
131 s - MPI_NUM_PROCS=12 AMBIENT_MKL_NUM_THREADS=2 run dmrg short 

74.4 s -                 AMBIENT_MKL_NUM_THREADS=1 run dmrg short 
82 s  -                  AMBIENT_MKL_NUM_THREADS=2 run dmrg short
86 s  -                  AMBIENT_MKL_NUM_THREADS=3 run dmrg short 
369.5 s -                                          run dmrg short # resource overloading


Ambient optional compiler defines:

AMBIENT_CHECK_BOUNDARIES     # checks memory boundaries safety every memptf calls
AMBIENT_NUMERIC_EXPERIMENTAL # enable experimental math-algorithms (i.e. svd or strassen gemm)
AMBIENT_LOOSE_FUTURE         # more aggressive lazy evaluation of futures (not safe)


Ambient concept caveats:

- The copy is performed by fusing two histories together (they share the same revision in one point). Therefore the element access for writing is unsafe (technically it will modify two objects at the same time).
