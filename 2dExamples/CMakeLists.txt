#ADD_SUBDIRECTORY(interface)
ADD_SUBDIRECTORY(dynamic)

# create symlinks to the necessary mesh files
execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_SOURCE_DIR}/myNutoExamples/meshFiles/2d/2d_miehe_symmetric_three_point_bending.msh" "${CMAKE_CURRENT_BINARY_DIR}/2d_miehe_symmetric_three_point_bending.msh"
)

# 2d_phase_field_3_point_miehe.cpp
ADD_EXECUTABLE(2d_phase_field_3_point_miehe 2d_phase_field_3_point_miehe.cpp)
TARGET_LINK_LIBRARIES(2d_phase_field_3_point_miehe NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
IF (ENABLE_VISUALIZE)
  TARGET_LINK_LIBRARIES(2d_phase_field_3_point_miehe NuToVisualize)
ENDIF (ENABLE_VISUALIZE)
IF(MUMPS_FOUND)
  TARGET_LINK_LIBRARIES(2d_phase_field_3_point_miehe ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
ENDIF(MUMPS_FOUND)

# 2d_local_damage_3_point_bending.cpp
ADD_EXECUTABLE(2d_local_damage_3_point_bending 2d_local_damage_3_point_bending.cpp)
TARGET_LINK_LIBRARIES(2d_local_damage_3_point_bending NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
IF (ENABLE_VISUALIZE)
  TARGET_LINK_LIBRARIES(2d_local_damage_3_point_bending NuToVisualize)
ENDIF (ENABLE_VISUALIZE)
IF(MUMPS_FOUND)
  TARGET_LINK_LIBRARIES(2d_local_damage_3_point_bending ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
ENDIF(MUMPS_FOUND)


## 2d_phase_field_single_edge_notched_tension.cpp
#ADD_EXECUTABLE(2d_phase_field_single_edge_notched_tension 2d_phase_field_single_edge_notched_tension.cpp)
#TARGET_LINK_LIBRARIES(2d_phase_field_single_edge_notched_tension NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
#IF (ENABLE_VISUALIZE)
#  TARGET_LINK_LIBRARIES(2d_phase_field_single_edge_notched_tension NuToVisualize)
#ENDIF (ENABLE_VISUALIZE)
#IF(MUMPS_FOUND)
#  TARGET_LINK_LIBRARIES(2d_phase_field_single_edge_notched_tension ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
#ENDIF(MUMPS_FOUND)

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

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})


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
  # 2d_feti_total_mpi_sparse_new_mesh_format
  ADD_EXECUTABLE(2d_feti_total_mpi_sparse_new_mesh_format 2d_feti_total_mpi_sparse_new_mesh_format.cpp)
  target_link_libraries(2d_feti_total_mpi_sparse_new_mesh_format ${MPI_LIBRARIES})
  TARGET_LINK_LIBRARIES(2d_feti_total_mpi_sparse_new_mesh_format NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(2d_feti_total_mpi_sparse_new_mesh_format NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(2d_feti_total_mpi_sparse_new_mesh_format ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(2d_feti_total_mpi_sparse_new_mesh_format PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)

