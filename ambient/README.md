Ambient
=======
**C++ framework for automatic application parallelisation under shared/distributed memory systems running Linux/OSX**

### Compilation defines

- AMBIENT_THREADING  
  set the desired threading implementation (CILK, OPENMP or SERIAL)  
  [auto]

- AMBIENT_MPI  
  MPI mode (use MPI_DISABLE or set the desired threading level)  
  [MPI_THREAD_FUNNELED]
                                            
- AMBIENT_DEFAULT_IB  
  Default blocking factor (partition/tile default size)  
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

### Environment variables

- AMBIENT_VERBOSE  
  print-out Ambient configuration prior to running  
  [not set]

- AMBIENT_BULK_LIMIT=[p]  
  limit the data bulk memory consumption by [p] percents of total memory  
  [60]
                                            
- AMBIENT_BULK_REUSE  
  setting this variable will enable bulk garbage collection  
  [not set]
                                            
- AMBIENT_BULK_FORCE_FREE  
  deallocate data bulk every time the sync has finished  
  [not set]

### License

Distributed under the Boost Software License, Version 1.0.  
(See http://www.boost.org/LICENSE_1_0.txt)
