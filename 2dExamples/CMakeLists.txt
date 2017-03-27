#ADD_SUBDIRECTORY(interface)

ADD_SUBDIRECTORY(feti)
ADD_SUBDIRECTORY(phase_field)
ADD_SUBDIRECTORY(benchmark)

set(MESH_DIR "${CMAKE_SOURCE_DIR}/applications/custom/myNutoExamples/meshFiles/2d")


execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/feti/fine/feti_beam.msh" "${CMAKE_CURRENT_BINARY_DIR}/feti_beam.msh"
)








  # 2d_feti_structure
  foreach(file
          #2d_feti_structure
          #2d_feti_structure_lu_solver
          )
    add_executable        (${file} ${file}.cpp)
    target_link_libraries (${file} ${MPI_LIBRARIES})
    target_link_libraries (${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    target_link_libraries (${file} NuToVisualize)
    target_link_libraries (${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  endforeach()

foreach(file
        #2d_feti_structure_compare_phase_field
        )
  add_executable        (${file} ${file}.cpp)
  target_link_libraries (${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
  target_link_libraries (${file} NuToVisualize)
  target_link_libraries (${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endforeach()

