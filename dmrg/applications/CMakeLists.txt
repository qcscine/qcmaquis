add_subdirectory(measure)

if (BUILD_INIT)
  add_subdirectory(init)
endif(BUILD_INIT)

if (BUILD_CIDEAS)
  add_subdirectory(cideas)
endif(BUILD_CIDEAS)

if (BUILD_MPS_ROTATE)
  add_subdirectory(rotate)
endif (BUILD_MPS_ROTATE)

if (BUILD_MPS_COMPRESS)
  # mps_compress does not work!!
  add_subdirectory(compress)
endif (BUILD_MPS_COMPRESS)

add_subdirectory(tools)

# Main executable
set(DMRG_APP_LIBRARIES maquis_dmrg dmrg_models dmrg_utils ${DMRG_LIBRARIES})
add_executable(qcmaquis qcmaquis.cpp)
target_link_libraries(qcmaquis ${DMRG_APP_LIBRARIES})
install(TARGETS qcmaquis RUNTIME DESTINATION bin COMPONENT applications)
