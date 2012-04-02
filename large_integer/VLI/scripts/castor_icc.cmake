set(CMAKE_C_COMPILER icc)
set(CMAKE_CXX_COMPILER icpc)
set (VLI_COMPILER_FLAGS_RELEASE " -Wall -O2 -funroll-loops -m64 -std=c++0x" CACHE STRING "Compiler flags for a regular compile.")
set (VLI_COMPILER_FLAGS_DEBUG " -Wall -g -O0 -m64" CACHE STRING "Compiler flags for a debug compile.")
