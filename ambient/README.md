Ambient
=======

### Installation
Ambient is a header-only library (so just installing the sources is enough):  
```sh
cmake .  
make install
```

### Usage
Compilation of the target application against Ambient's include folder with C++11 enabled is generally sufficient. Otherwise be sure to check the compilation options: threading backend and MPI mode are compiler-specific by default (so the respective compiler knobs might be needed). 

To enforce threading backend or MPI mode use the following knobs:

    -DAMBIENT_THREADING [CILK, OPENMP or NONE]  
    -DAMBIENT_MPI [MPI_DISABLE or the desired MPI threading level]

(Set AMBIENT_VERBOSE environment variable to see the resulting configuration).

Check the [developer's guide](http://ambient.comp-phys.org/guide.pdf) for more detailed information.

### License
    Distributed under the Boost Software License, Version 1.0.  
    (See http://www.boost.org/LICENSE_1_0.txt)
