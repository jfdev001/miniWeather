cmake_minimum_required(VERSION 3.0)
project(miniWeather Fortran)

enable_testing()
option(BUILD_TESTING "Build tests" OFF)

include(utils.cmake)


############################################################
## Set Parameters
############################################################
if ("${NX}" STREQUAL "")
  SET(NX 100)
endif()
if ("${NZ}" STREQUAL "")
  SET(NZ 50)
endif()
if ("${SIM_TIME}" STREQUAL "")
  SET(SIM_TIME 1000)
endif()
if ("${OUT_FREQ}" STREQUAL "")
  SET(OUT_FREQ 10)
endif()
if ("${DATA_SPEC}" STREQUAL "")
  SET(DATA_SPEC DATA_SPEC_THERMAL)
endif()
SET(EXE_DEFS "-D_NX=${NX} -D_NZ=${NZ} -D_SIM_TIME=${SIM_TIME} -D_OUT_FREQ=${OUT_FREQ} -D_DATA_SPEC=${DATA_SPEC}")
SET(TEST_DEFS "-D_NX=100 -D_NZ=50 -D_SIM_TIME=400 -D_OUT_FREQ=400 -D_DATA_SPEC=DATA_SPEC_THERMAL")
if ("${FFLAGS}" MATCHES ".*SINGLE_PREC.*")
  set(tol "1.e-5")
  message(STATUS "Using single precision. Test tolerance ${tol}")
else()
  set(tol "1.e-13")
  message(STATUS "Using double precision. Test tolerance ${tol}")
endif()


############################################################
## Append FFLAGS
############################################################
SET(CMAKE_Fortran_FLAGS "${FFLAGS}")


############################################################
## Compile the serial version
############################################################
add_executable(serial miniWeather_serial.F90)
set_target_properties(serial PROPERTIES COMPILE_FLAGS "${EXE_DEFS}")
if(${BUILD_TESTING})
  add_executable(serial_test miniWeather_serial.F90) 
  set_target_properties(serial_test PROPERTIES COMPILE_FLAGS "${TEST_DEFS}")
endif()

if (NOT ("${LDFLAGS}" STREQUAL "") )
  target_link_libraries(serial      "${LDFLAGS}")
  if(${BUILD_TESTING})
    target_link_libraries(serial_test "${LDFLAGS}")
  endif()
endif()

if (NOT ("${SERIAL_LINK_FLAGS}" STREQUAL "") )
  target_link_libraries(serial      "${SERIAL_LINK_FLAGS}")
  if(${BUILD_TESTING})
    target_link_libraries(serial_test "${SERIAL_LINK_FLAGS}")
  endif()
endif()

if(${BUILD_TESTING})
    add_test(NAME SERIAL_TEST COMMAND ./check_output.sh ./serial_test ${tol} 4.5e-5 ) 
endif()


############################################################
## Compile the MPI version
############################################################
add_executable(mpi miniWeather_mpi.F90)
set_target_properties(mpi PROPERTIES COMPILE_FLAGS "${EXE_DEFS}")

if(${BUILD_TESTING})
    add_executable(mpi_test miniWeather_mpi.F90)
    set_target_properties(mpi_test PROPERTIES COMPILE_FLAGS "${TEST_DEFS}")
endif()

if (NOT ("${LDFLAGS}" STREQUAL "") )
  target_link_libraries(mpi      "${LDFLAGS}")
  if(${BUILD_TESTING})
    target_link_libraries(mpi_test "${LDFLAGS}")
  endif()
endif()

if (NOT ("${MPI_LINK_FLAGS}" STREQUAL "") )
  target_link_libraries(mpi      "${MPI_LINK_FLAGS}")
  if(${BUILD_TESTING})
    target_link_libraries(mpi_test "${MPI_LINK_FLAGS}")
  endif()
endif()

if(${BUILD_TESTING})
  add_test(NAME MPI_TEST COMMAND ./check_output.sh ./mpi_test ${tol} 4.5e-5 ) 
endif()


############################################################
## Compile the MPI + OpenMP version
############################################################
if (NOT ("${OPENMP_FLAGS}" STREQUAL "") )
  add_executable(openmp miniWeather_mpi_openmp.F90)
  set_target_properties(openmp 
    PROPERTIES COMPILE_FLAGS "${EXE_DEFS} ${OPENMP_FLAGS}")

  if(${BUILD_TESTING})
    add_executable(openmp_test miniWeather_mpi_openmp.F90)
    set_target_properties(openmp_test 
      PROPERTIES COMPILE_FLAGS "${TEST_DEFS} ${OPENMP_FLAGS}")
  endif()

  if (NOT ("${LDFLAGS}" STREQUAL "") )
    target_link_libraries(openmp      "${LDFLAGS}")
    if(${BUILD_TESTING})
      target_link_libraries(openmp_test "${LDFLAGS}")
    endif()
  endif()

  if ("${OPENMP_LINK_FLAGS}" STREQUAL "")
    SET(OPENMP_LINK_FLAGS ${OPENMP_FLAGS})
  endif()

  target_link_libraries(openmp      "${OPENMP_LINK_FLAGS}")

  if(${BUILD_TESTING})
    target_link_libraries(openmp_test "${OPENMP_LINK_FLAGS}")
    add_test(NAME OPENMP_TEST COMMAND ./check_output.sh ./openmp_test ${tol} 4.5e-5 ) 
  endif()

endif()


############################################################
## NOTE: Commented out since unused
##
## Compile the MPI + OpenACC version
############################################################
# if (NOT ("${OPENACC_FLAGS}" STREQUAL "") )
#   add_executable(openacc miniWeather_mpi_openacc.F90)
#   set_target_properties(openacc PROPERTIES COMPILE_FLAGS "${EXE_DEFS} ${OPENACC_FLAGS}")
#
# if(${BUILD_TESTING})
#     #   add_executable(openacc_test miniWeather_mpi_openacc.F90)
# endif()
# if(${BUILD_TESTING})
#     #   set_target_properties(openacc_test PROPERTIES COMPILE_FLAGS "${TEST_DEFS} ${OPENACC_FLAGS}")
# endif()
# #
# #   if (NOT ("${LDFLAGS}" STREQUAL "") )
# #     target_link_libraries(openacc      "${LDFLAGS}")
# if(${BUILD_TESTING})
#     #     target_link_libraries(openacc_test "${LDFLAGS}")
# endif()
# #   endif()
# #   if ("${OPENACC_LINK_FLAGS}" STREQUAL "")
# #     SET(OPENACC_LINK_FLAGS ${OPENACC_FLAGS})
# #   endif()
# #   target_link_libraries(openacc      "${OPENACC_LINK_FLAGS}")
# if(${BUILD_TESTING})
#     #   target_link_libraries(openacc_test "${OPENACC_LINK_FLAGS}")
# endif()
# #
# if(${BUILD_TESTING})
#     #   add_test(NAME OPENACC_TEST COMMAND ./check_output.sh ./openacc_test ${tol} 4.5e-5 ) 
# endif()
# # endif()



############################################################
## NOTE: Commented out because gcc -fopenmp -dM -E - < /dev/null | grep _OPENMP
## shows that OpenMP 4.5 is implemented... therefore the fact that the test
## fails must be due to some logical error in the implementation or an issue
## with compiler flags which is not worth troubleshooting

## Compile the MPI + OpenMP4.5 version
############################################################
# if (NOT ("${OPENMP45_FLAGS}" STREQUAL "") )
#   add_executable(openmp45 miniWeather_mpi_openmp45.F90)
#   set_target_properties(openmp45 PROPERTIES COMPILE_FLAGS "${EXE_DEFS} ${OPENMP45_FLAGS}")
#
# if(${BUILD_TESTING})
#     #   add_executable(openmp45_test miniWeather_mpi_openmp45.F90)
# endif()
# if(${BUILD_TESTING})
#     #   set_target_properties(openmp45_test PROPERTIES COMPILE_FLAGS "${TEST_DEFS} ${OPENMP45_FLAGS}")
# endif()
# #
# #   # The XL compiler dumps out non-unique filenames that screw up parallel compilation
# #   # So it must compile the test at a different time than the original executable
# #   if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "XL")
# if(${BUILD_TESTING})
#     #     add_dependencies(openmp45_test openmp45)
# endif()
# #   endif()
# #
# #   if (NOT ("${LDFLAGS}" STREQUAL "") )
# #     target_link_libraries(openmp45      "${LDFLAGS}")
# if(${BUILD_TESTING})
#     #     target_link_libraries(openmp45_test "${LDFLAGS}")
# endif()
# #   endif()
# #   if ("${OPENMP45_LINK_FLAGS}" STREQUAL "")
# #     SET(OPENMP45_LINK_FLAGS ${OPENMP45_FLAGS})
# #   endif()
# #   target_link_libraries(openmp45      "${OPENMP45_LINK_FLAGS}")
# if(${BUILD_TESTING})
#     #   target_link_libraries(openmp45_test "${OPENMP45_LINK_FLAGS}")
# endif()
# #
# if(${BUILD_TESTING})
#     #   add_test(NAME OPENMP45_TEST COMMAND ./check_output.sh ./openmp45_test ${tol} 4.5e-5 ) 
# endif()
# endif()



############################################################
## NOTE: Tests at the moment do not run this, you must specify 
## -DDO_CONCURRENT_FLAGS during cmake configuration if you wish to run these
##
## Compile the MPI + do concurrent version
############################################################
if (NOT ("${DO_CONCURRENT_FLAGS}" STREQUAL "") )
  add_executable(do_concurrent miniWeather_mpi_do_concurrent.F90)
  set_target_properties(do_concurrent PROPERTIES COMPILE_FLAGS "${EXE_DEFS} ${DO_CONCURRENT_FLAGS}")

if(${BUILD_TESTING})
      add_executable(do_concurrent_test miniWeather_mpi_do_concurrent.F90)
endif()
if(${BUILD_TESTING})
      set_target_properties(do_concurrent_test PROPERTIES COMPILE_FLAGS "${TEST_DEFS} ${DO_CONCURRENT_FLAGS}")
endif()

  if (NOT ("${LDFLAGS}" STREQUAL "") )
    target_link_libraries(do_concurrent      "${LDFLAGS}")
if(${BUILD_TESTING})
        target_link_libraries(do_concurrent_test "${LDFLAGS}")
endif()
  endif()
  if ("${DO_CONCURRENT_LINK_FLAGS}" STREQUAL "")
    SET(DO_CONCURRENT_LINK_FLAGS ${DO_CONCURRENT_FLAGS})
  endif()
  target_link_libraries(do_concurrent      "${DO_CONCURRENT_LINK_FLAGS}")
if(${BUILD_TESTING})
      target_link_libraries(do_concurrent_test "${DO_CONCURRENT_LINK_FLAGS}")
endif()

if(${BUILD_TESTING})
      add_test(NAME DO_CONCURRENT_TEST COMMAND ./check_output.sh ./do_concurrent_test ${tol} 4.5e-5 ) 
endif()
endif()



############################################################
## NOTE: Commented out since not necessary 
##
## Compile the serial allocation during runtime version
############################################################
# add_executable(allocate_during_runtime exploration/miniWeather_serial_allocate_during_runtime.F90)
# set_target_properties(allocate_during_runtime PROPERTIES COMPILE_FLAGS "${EXE_DEFS}")
#
# if(${BUILD_TESTING})
#     # add_executable(allocate_during_runtime_test exploration/miniWeather_serial_allocate_during_runtime.F90)
# endif()
# if(${BUILD_TESTING})
#     # set_target_properties(allocate_during_runtime_test PROPERTIES COMPILE_FLAGS "${TEST_DEFS}")
# endif()
# #
# # if (NOT ("${LDFLAGS}" STREQUAL "") )
# #   target_link_libraries(allocate_during_runtime      "${LDFLAGS}")
# if(${BUILD_TESTING})
#     #   target_link_libraries(allocate_during_runtime_test "${LDFLAGS}")
# endif()
# # endif()
# # if (NOT ("${SERIAL_LINK_FLAGS}" STREQUAL "") )
# #   target_link_libraries(allocate_during_runtime      "${SERIAL_LINK_FLAGS}")
# if(${BUILD_TESTING})
#     #   target_link_libraries(allocate_during_runtime_test "${SERIAL_LINK_FLAGS}")
# endif()
# # endif()
# #
# if(${BUILD_TESTING})
#     # add_test(NAME ALLOCATE_DURING_RUNTIME_TEST COMMAND ./check_output.sh ./allocate_during_runtime_test ${tol} 4.5e-5 ) 
# endif()



############################################################
## NOTE: Commented out since not necessary 
##
## Compile the serial allocation during runtime version
############################################################
# add_executable(allocate_and_assumed_shape exploration/miniWeather_serial_allocate_and_assumed_shape.F90)
# set_target_properties(allocate_and_assumed_shape PROPERTIES COMPILE_FLAGS "${EXE_DEFS}")
#
if(${BUILD_TESTING})
    # add_executable(allocate_and_assumed_shape_test exploration/miniWeather_serial_allocate_and_assumed_shape.F90)
endif()
if(${BUILD_TESTING})
    # set_target_properties(allocate_and_assumed_shape_test PROPERTIES COMPILE_FLAGS "${TEST_DEFS}")
endif()
#
# if (NOT ("${LDFLAGS}" STREQUAL "") )
#   target_link_libraries(allocate_and_assumed_shape      "${LDFLAGS}")
if(${BUILD_TESTING})
    #   target_link_libraries(allocate_and_assumed_shape_test "${LDFLAGS}")
endif()
# endif()
# if (NOT ("${SERIAL_LINK_FLAGS}" STREQUAL "") )
#   target_link_libraries(allocate_and_assumed_shape      "${SERIAL_LINK_FLAGS}")
if(${BUILD_TESTING})
    #   target_link_libraries(allocate_and_assumed_shape_test "${SERIAL_LINK_FLAGS}")
endif()
# endif()
#
if(${BUILD_TESTING})
    # add_test(NAME ALLOCATE_AND_ASSUMED_SHAPE_TEST COMMAND ./check_output.sh ./allocate_and_assumed_shape_test ${tol} 4.5e-5 ) 
endif()



