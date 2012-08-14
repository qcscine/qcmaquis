#set(COMPILER_FLAGS_VAR CMAKE_CXX_FLAGS_RELEASE) # Quite just the release mode
set(COMPILER_FLAGS_VAR CMAKE_CXX_FLAGS) # Quite everything

set(${COMPILER_FLAGS_VAR} "-Wno-return-type-c-linkage ${${COMPILER_FLAGS_VAR}}")
set(${COMPILER_FLAGS_VAR} "-Wno-dangling-else ${${COMPILER_FLAGS_VAR}}")
