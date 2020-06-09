# ALPS CMakeLists.txt modified for distribution with QCMaquis
# (c) 2020, Leon Freitag
# Original copyright header:
#
#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#
# libalps
#

set(ALPS_SOURCES "")
set(ALPS_SHARED_CPPFLAGS ALPS_EXPORTS=1)
set(ALPS_STATIC_CPPFLAGS "")
set(BUILD_SHARED_LIBS OFF)
set(ALPS_SOURCES
    # utility
    utility/encode.cpp utility/copyright.cpp utility/os.cpp utility/temporary_filename.cpp
    utility/vmusage.cpp
#     # alea
#     alea/nan.C alea/observable.C alea/observablefactory.C alea/observableset.C
    # expression
    expression/evaluator.C
    # # lattice
    # lattice/disorder.C lattice/latticedescriptor.C lattice/latticegraphdescriptor.C
    # lattice/latticelibrary.C lattice/unitcell.C
    # # model
    # model/bondterm.C model/globaloperator.C model/modellibrary.C model/siteterm.C
    # osiris
    osiris/buffer.C osiris/dump.C osiris/dumparchive.C osiris/process.C osiris/comm.C
    osiris/mpdump.C osiris/xdrdump.C osiris/xdrcore.C
    # parameter
    parameter/parameter.C parameter/parameterlist.C parameter/parameters.C
    # parser
    parser/parser.C parser/xmlattributes.C parser/xmlhandler.C parser/xmlparser.C
    parser/xmlstream.C parser/xslt_path.C
    #hdf5
    hdf5/archive.cpp ngs/lib/signal.cpp ngs/lib/stacktrace.cpp
    # random
    random/rngfactory.C random/parallel/get_prime.cpp random/parallel/primelist_64.cpp
#     # scheduler
#     scheduler/abstract_task.C scheduler/convertxdr.C scheduler/factory.C scheduler/info.C
#     scheduler/master_scheduler.C scheduler/measurement_operators.C
#     scheduler/montecarlo.C scheduler/mpp_scheduler.C scheduler/options.C
#     scheduler/remote_task.C scheduler/remote_worker.C scheduler/scheduler.C
#     scheduler/serial_scheduler.C scheduler/signal.C scheduler/single_scheduler.C
#     scheduler/slave_task.C scheduler/task.C scheduler/worker.C scheduler/workertask.C
#     # ngs sources
#     stop_callback.cpp ngs/lib/api.cpp ngs/lib/short_print.cpp ngs/lib/sleep.cpp
    # `
    # parameters
    # ngs/lib/paramvalue.cpp ngs/lib/paramproxy.cpp ngs/lib/params.cpp
#     # scheduler
#     mcbase.cpp parseargs.cpp
# #deprecated, TODO: replace by new alea
#     ngs/lib/mcoptions.cpp ngs/lib/observablewrappers.cpp ngs/lib/make_deprecated_parameters.cpp ngs/lib/make_parameters_from_xml.cpp
#     ngs/lib/mcresult.cpp ngs/lib/mcresults.cpp ngs/lib/mcobservables.cpp ngs/lib/mcobservable.cpp
  )

#   if (ALPS_BUILD_PYTHON)
#     set(ALPS_SOURCES ${ALPS_SOURCES} ngs/lib/get_numpy_type.cpp hdf5/python.cpp python/numpy_array.cpp)
#   endif (ALPS_BUILD_PYTHON)

#   # OpenMPI ULFM
#   if (ALPS_NGS_OPENMPI_ULFM)
#     set(ALPS_SOURCES ${ALPS_SOURCES} ngs/lib/ulfm.cpp)
#   endif (ALPS_NGS_OPENMPI_ULFM)

#   # new alea
#   if (ALPS_NGS_USE_NEW_ALEA)
# #    set(ALPS_SOURCES ${ALPS_SOURCES} ngs/lib/alea/accumulator_set.cpp ngs/lib/alea/result_set.cpp)
#     set(ALPS_SOURCES ${ALPS_SOURCES} ngs/accumulator/accumulator.cpp)
#   endif (ALPS_NGS_USE_NEW_ALEA)

# if (NOT ALPS_FOR_VISTRAILS)
#   set(ALPS_SOURCES ${ALPS_SOURCES}
#       # parapack
#       parapack/clone.C parapack/clone_info.C parapack/filelock.C parapack/logger.C
#       parapack/job.C parapack/mc_worker.C parapack/measurement.C parapack/process_impl.C
#       parapack/option.C parapack/parapack.C parapack/queue.C parapack/rng_helper.C
#       parapack/types.C parapack/util.C parapack/version.C parapack/worker_factory.C
#       parapack/staging.C
#       # ngs-parapack
#       ngs/lib/clone.cpp ngs/lib/clone_info.cpp ngs/lib/job.cpp ngs/lib/parapack.cpp ngs/lib/worker_factory.cpp
#   )
# endif(NOT ALPS_FOR_VISTRAILS)

add_library(alps SHARED ${ALPS_SOURCES})
set_property(TARGET alps PROPERTY POSITION_INDEPENDENT_CODE TRUE)

install(TARGETS alps EXPORT DMRGTargets COMPONENT libraries DESTINATION lib)
# set(ALPS_HEADER_DIR "include")
#     install(DIRECTORY src/alps src/boost src/ietl src/mocasito COMPONENT headers DESTINATION ${ALPS_HEADER_DIR}
#             FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.ipp"
#             PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
#     install(FILES ${PROJECT_BINARY_DIR}/src/alps/config.h ${PROJECT_BINARY_DIR}/src/alps/version.h DESTINATION ${ALPS_HEADER_DIR}/alps COMPONENT headers)
export(TARGETS alps APPEND FILE "${PROJECT_BINARY_DIR}/DMRGTargets.cmake")
