set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set (VLI_COMPILER_FLAGS_RELEASE " -Wall -Os -funroll-loops -m64" CACHE STRING "Compiler flags for a regular compile.")
set (VLI_COMPILER_FLAGS_DEBUG " -Wall -g -O0 -m64" CACHE STRING "Compiler flags for a debug compile.")
