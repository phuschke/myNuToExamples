#ADD_SUBDIRECTORY(interface)
ADD_SUBDIRECTORY(dynamic)

# create symlinks to the necessary mesh files
execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_SOURCE_DIR}/myNutoExamples/meshFiles/2d/2d_miehe_symmetric_three_point_bending.msh" "${CMAKE_CURRENT_BINARY_DIR}/2d_miehe_symmetric_three_point_bending.msh"
)

# create symlinks to the necessary mesh files
execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_SOURCE_DIR}/myNutoExamples/meshFiles/2d/feti/fine/feti.msh_000001" "${CMAKE_CURRENT_BINARY_DIR}/feti.msh_000001"
)
execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_SOURCE_DIR}/myNutoExamples/meshFiles/2d/feti/fine/feti.msh_000002" "${CMAKE_CURRENT_BINARY_DIR}/feti.msh_000002"
)
execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_SOURCE_DIR}/myNutoExamples/meshFiles/2d/feti/fine/feti.msh_000003" "${CMAKE_CURRENT_BINARY_DIR}/feti.msh_000003"
)

execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_SOURCE_DIR}/myNutoExamples/meshFiles/2d/feti/fine/feti_beam.msh" "${CMAKE_CURRENT_BINARY_DIR}/feti_beam.msh"
)

# read_json.cpp
ADD_EXECUTABLE(read_json read_json.cpp jsoncpp.cpp)
TARGET_LINK_LIBRARIES(read_json NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
IF (ENABLE_VISUALIZE)
  TARGET_LINK_LIBRARIES(read_json NuToVisualize)
ENDIF (ENABLE_VISUALIZE)
IF(MUMPS_FOUND)
  TARGET_LINK_LIBRARIES(read_json ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
ENDIF(MUMPS_FOUND)

# 2d_feti_structure_compare_phase_field.cpp
ADD_EXECUTABLE(2d_feti_structure_compare_phase_field 2d_feti_structure_compare_phase_field.cpp)
TARGET_LINK_LIBRARIES(2d_feti_structure_compare_phase_field NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
IF (ENABLE_VISUALIZE)
  TARGET_LINK_LIBRARIES(2d_feti_structure_compare_phase_field NuToVisualize)
ENDIF (ENABLE_VISUALIZE)
IF(MUMPS_FOUND)
  TARGET_LINK_LIBRARIES(2d_feti_structure_compare_phase_field ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
ENDIF(MUMPS_FOUND)

# 2d_feti_structure_compare_gradient_damage.cpp
ADD_EXECUTABLE(2d_feti_structure_compare_gradient_damage 2d_feti_structure_compare_gradient_damage.cpp)
TARGET_LINK_LIBRARIES(2d_feti_structure_compare_gradient_damage NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
IF (ENABLE_VISUALIZE)
  TARGET_LINK_LIBRARIES(2d_feti_structure_compare_gradient_damage NuToVisualize)
ENDIF (ENABLE_VISUALIZE)
IF(MUMPS_FOUND)
  TARGET_LINK_LIBRARIES(2d_feti_structure_compare_gradient_damage ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
ENDIF(MUMPS_FOUND)



# 2d_feti_structure_compare.cpp
#ADD_EXECUTABLE(2d_feti_structure_compare 2d_feti_structure_compare.cpp)
#TARGET_LINK_LIBRARIES(2d_feti_structure_compare NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
#IF (ENABLE_VISUALIZE)
#  TARGET_LINK_LIBRARIES(2d_feti_structure_compare NuToVisualize)
#ENDIF (ENABLE_VISUALIZE)
#IF(MUMPS_FOUND)
#  TARGET_LINK_LIBRARIES(2d_feti_structure_compare ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
#ENDIF(MUMPS_FOUND)

# 2d_local_damage_3_point_bending.cpp
#ADD_EXECUTABLE(2d_local_damage_3_point_bending 2d_local_damage_3_point_bending.cpp)
#TARGET_LINK_LIBRARIES(2d_local_damage_3_point_bending NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
#IF (ENABLE_VISUALIZE)
#  TARGET_LINK_LIBRARIES(2d_local_damage_3_point_bending NuToVisualize)
#ENDIF (ENABLE_VISUALIZE)
#IF(MUMPS_FOUND)
#  TARGET_LINK_LIBRARIES(2d_local_damage_3_point_bending ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
#ENDIF(MUMPS_FOUND)


# 2d_phase_field_single_edge_notched_tension.cpp
ADD_EXECUTABLE(2d_phase_field_single_edge_notched_tension 2d_phase_field_single_edge_notched_tension.cpp)
TARGET_LINK_LIBRARIES(2d_phase_field_single_edge_notched_tension NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
IF (ENABLE_VISUALIZE)
  TARGET_LINK_LIBRARIES(2d_phase_field_single_edge_notched_tension NuToVisualize)
ENDIF (ENABLE_VISUALIZE)
IF(MUMPS_FOUND)
  TARGET_LINK_LIBRARIES(2d_phase_field_single_edge_notched_tension ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
ENDIF(MUMPS_FOUND)

## 2d_phase_field_single_edge_notched_shear.cpp
#ADD_EXECUTABLE(2d_phase_field_single_edge_notched_shear 2d_phase_field_single_edge_notched_shear.cpp)
#TARGET_LINK_LIBRARIES(2d_phase_field_single_edge_notched_shear NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
#IF (ENABLE_VISUALIZE)
#TARGET_LINK_LIBRARIES(2d_phase_field_single_edge_notched_shear NuToVisualize)
#ENDIF (ENABLE_VISUALIZE)
#IF(MUMPS_FOUND)
#TARGET_LINK_LIBRARIES(2d_phase_field_single_edge_notched_shear ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
#ENDIF(MUMPS_FOUND)


message(status "##########################################")

#find_package(MPI REQUIRED)
#include_directories(${MPI_INCLUDE_PATH})


#if(MPI_COMPILE_FLAGS)
#  set_target_properties(hello PROPERTIES
#    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
#endif()
#
#if(MPI_LINK_FLAGS)
#  set_target_properties(hello PROPERTIES
#    LINK_FLAGS "${MPI_LINK_FLAGS}")
#endif()
message("{MPI_INCLUDE_PATH}" ${MPI_INCLUDE_PATH})
message("{MPI_COMPILE_FLAGS}" ${MPI_COMPILE_FLAGS})
message("{MPI_LINK_FLAGS}" ${MPI_LINK_FLAGS})
message("{MPI_LIBRARIES}" ${MPI_LIBRARIES})

message(status "##########################################")
  # 2d_feti_total_mpi_sparse
 # ADD_EXECUTABLE(2d_feti_total_mpi_sparse 2d_feti_total_mpi_sparse.cpp)
 # target_link_libraries(2d_feti_total_mpi_sparse ${MPI_LIBRARIES})
 # TARGET_LINK_LIBRARIES(2d_feti_total_mpi_sparse NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
 # IF (ENABLE_VISUALIZE)
 #   TARGET_LINK_LIBRARIES(2d_feti_total_mpi_sparse NuToVisualize)
 # ENDIF (ENABLE_VISUALIZE)
 # IF(MUMPS_FOUND)
 #   TARGET_LINK_LIBRARIES(2d_feti_total_mpi_sparse ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
 # ENDIF(MUMPS_FOUND)
 # IF(MINGW)
 #   SET_TARGET_PROPERTIES(2d_feti_total_mpi_sparse PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
 # ENDIF(MINGW)


  # 2d_test
 # ADD_EXECUTABLE(2d_test 2d_test.cpp)
 # target_link_libraries(2d_test ${MPI_LIBRARIES})
 # TARGET_LINK_LIBRARIES(2d_test NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
 # IF (ENABLE_VISUALIZE)
 #   TARGET_LINK_LIBRARIES(2d_test NuToVisualize)
 # ENDIF (ENABLE_VISUALIZE)
 # IF(MUMPS_FOUND)
 #   TARGET_LINK_LIBRARIES(2d_test ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
 # ENDIF(MUMPS_FOUND)
 # IF(MINGW)
 #   SET_TARGET_PROPERTIES(2d_test PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
 # ENDIF(MINGW)
 #
 #
  # 2d_feti_structure_phase_field
  ADD_EXECUTABLE(2d_feti_structure_phase_field 2d_feti_structure_phase_field.cpp)
  target_link_libraries(2d_feti_structure_phase_field ${MPI_LIBRARIES})
  TARGET_LINK_LIBRARIES(2d_feti_structure_phase_field NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(2d_feti_structure_phase_field NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(2d_feti_structure_phase_field ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(2d_feti_structure_phase_field PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)


  # 2d_feti_structure
  ADD_EXECUTABLE(2d_feti_structure 2d_feti_structure.cpp)
  target_link_libraries(2d_feti_structure ${MPI_LIBRARIES})
  TARGET_LINK_LIBRARIES(2d_feti_structure NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(2d_feti_structure NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(2d_feti_structure ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(2d_feti_structure PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)

  # 2d_feti_structure_gradient_damage
  ADD_EXECUTABLE(2d_feti_structure_gradient_damage 2d_feti_structure_gradient_damage.cpp)
  target_link_libraries(2d_feti_structure_gradient_damage ${MPI_LIBRARIES})
  TARGET_LINK_LIBRARIES(2d_feti_structure_gradient_damage NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(2d_feti_structure_gradient_damage NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(2d_feti_structure_gradient_damage ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(2d_feti_structure_gradient_damage PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)
