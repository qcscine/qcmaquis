Ambient
=======
**Dataflow C++ framework for distributed computations**

### Ambient environment variables

- AMBIENT_VERBOSE  
  print-out Ambient configuration prior to running  
  [not set]

- AMBIENT_MKL_NUM_THREADS=[thread count]  
  enable selective threading: MKL will use only 1 thread unless sync is called passing mkl_parallel() - then the [thread count] is used.  
  [not set]

- MKL_NUM_THREADS=[thread count]  
  MKL will use [thread count] in all of its calls.  
  *Warning: can cause performance degradation due to resources overloading. Use AMBIENT_MKL_NUM_THREADS to override it.*  
  [auto]
                                            
- AMBIENT_BULK_LIMIT=[p]  
  limit the data bulk memory consumption by [p] percents of total memory  
  [60]
                                            
- AMBIENT_BULK_REUSE  
  setting this variable will enable bulk garbage collection  
  [not set]
                                            
- AMBIENT_FORCE_BULK_DEALLOCATION  
  deallocate data bulk every time the sync has finished  
  [not set]


### Ambient compilation defines

- [AMBIENT_CILK, AMBIENT_OMP, AMBIENT_SERIAL]  
  manually set the desired threading implementation  
  [compiler dependent]

- AMBIENT_MPI_THREADING  
  desired level of MPI threading (note: Ambient calls MPI routines through the main thread)  
  [MPI_THREAD_FUNNELED]
                                            
- AMBIENT_IB  
  Ambient internal block-size (also used as a matrix tile size)  
  [2048]
                                            
- AMBIENT_INSTR_BULK_CHUNK  
  size (bytes) of memory chunks for operations logging (async calls info)  
  [16MB]
                                            
- AMBIENT_DATA_BULK_CHUNK  
  size (bytes) of memory chunks for communications and temporary objects (> size of corresponding tiles)  
  [64MB]
                                            
- AMBIENT_SERIAL_COLLECTION  
  enable to make operations collection not thread-safe
  [not set]
                                            
- AMBIENT_MEMPTF_CHECK_BOUNDARIES  
  checks memory boundaries overflow in every memptf call (used for 2D memory copies)  
  [not set]


### Ambient implementation caveats

- *Direct element access is slow and should be used only for debugging.*

- *The copy is performed by fusing two version-histories together (they share the same revision in one point).
  Therefore the direct element access (note: deprecated, see above) for writing is unsafe.*
