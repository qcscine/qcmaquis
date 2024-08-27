#
# Temporary Copy paste from scine_dev-utils
#
macro(import_pybind11)
  # If the target already exists, do nothing
  if(NOT TARGET pybind11::pybind11)
    find_package(pybind11 2.10.4 EXACT QUIET)
    if(TARGET pybind11::pybind11)
      message(STATUS "Found pybind11 at ${pybind11_DIR}")
    else()
      # Download it instead
      include(DownloadProject)
      download_project(
        PROJ pybind11
        GIT_REPOSITORY      https://github.com/pybind/pybind11.git
        GIT_TAG             v2.10.4
        QUIET
      )
      add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
      # Final check if all went well
      if(EXISTS "${pybind11_SOURCE_DIR}/CMakeLists.txt")
        message(STATUS
          "Pybind11 was not found in your PATH, so it was downloaded."
        )
      else()
        string(CONCAT error_msg
          "Pybind11 could not be downloaded."
        )
        message(FATAL_ERROR ${error_msg})
      endif()
    endif()
  endif()
endmacro()
