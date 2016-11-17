#ADD_SUBDIRECTORY(interface)

ADD_SUBDIRECTORY(feti)

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





message(status "##########################################")

message("{MPI_INCLUDE_PATH}" ${MPI_INCLUDE_PATH})
message("{MPI_COMPILE_FLAGS}" ${MPI_COMPILE_FLAGS})
message("{MPI_LINK_FLAGS}" ${MPI_LINK_FLAGS})
message("{MPI_LIBRARIES}" ${MPI_LIBRARIES})

message(status "##########################################")





  # create symlinks to the necessary mesh files
  foreach(task RANGE 0 6)
      message("Task ${task}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "/home/phuschke/gmshReader_project/build/examples/feti_beam_fine.mesh${task}" "${CMAKE_CURRENT_BINARY_DIR}/feti_json.msh_${task}"
  )
    endforeach()


  # 2d_feti_structure
  set(file 2d_feti_structure)
  add_executable        (${file} ${file}.cpp)
  target_link_libraries (${file} ${MPI_LIBRARIES})
  target_link_libraries (${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
  target_link_libraries (${file} NuToVisualize)
  target_link_libraries (${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

  # 2d_feti_structure_gradient_damage
  set(file 2d_feti_structure_gradient_damage)
  add_executable        (${file} ${file}.cpp)
  target_link_libraries (${file} ${MPI_LIBRARIES})
  target_link_libraries (${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
  target_link_libraries (${file} NuToVisualize)
  target_link_libraries (${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
