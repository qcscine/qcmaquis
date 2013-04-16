
configure_file(types/version.h.in ${CMAKE_CURRENT_BINARY_DIR}/types/version.h)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/types/version.h
        DESTINATION include/types COMPONENT headers)

install(DIRECTORY types utils
        DESTINATION include
        COMPONENT headers
        ${EXCLUDE_FROM_MATRICES}
        PATTERN ".svn" EXCLUDE
        PATTERN ".DS_Store" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.C" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.cu" EXCLUDE
        PATTERN "*.in" EXCLUDE
        PATTERN "CMakeLists.txt" EXCLUDE
        PATTERN "Makefile*" EXCLUDE
	)

