set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
#acml serial lacml, multi proc lacml_mp          
set(DEFAULT_BLAS_LAPACK manual)
set(BLAS_LAPACK_MANUAL_LIBS -lacml_mv -lacml)
set(BLAS_LAPACK_MANUAL_LIBS_DIR /opt/acml/4.4.0/gfortran64/lib) 
set(BLAS_LAPACK_MANUAL_INCLUDES /opt/acml/4.4.0/gfortran64/include)
