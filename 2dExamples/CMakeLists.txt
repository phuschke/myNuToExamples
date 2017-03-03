#ADD_SUBDIRECTORY(interface)

ADD_SUBDIRECTORY(feti)
ADD_SUBDIRECTORY(phase_field)
ADD_SUBDIRECTORY(benchmark)

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








  # create symlinks to the necessary mesh files
  foreach(task RANGE 0 1)
      message("Task ${task}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "/home/phuschke/gmshReader_project/build/examples/feti_beam_fine_2_subdomains.mesh${task}" "${CMAKE_CURRENT_BINARY_DIR}/feti_beam_fine_2_subdomains.mesh${task}"
  )
    endforeach()


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
