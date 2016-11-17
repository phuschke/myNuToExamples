


# create symlinks to the necessary mesh files
foreach(task RANGE 0 4)
    message("Task ${task}")
execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "/home/phuschke/gmshReader_project/build/examples/fourPointBending.mesh${task}" "${CMAKE_CURRENT_BINARY_DIR}/fourPointBending.msh_${task}"
)
endforeach()

foreach(task RANGE 0 7)
    message("Task ${task}")
execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "/home/phuschke/gmshReader_project/build/examples/threePointBending.mesh${task}" "${CMAKE_CURRENT_BINARY_DIR}/threePointBending.msh_${task}"
)
endforeach()


foreach(task RANGE 0 4)
    message("Task ${task}")
execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "/home/phuschke/gmshReader_project/build/examples/feti_beam_fine_tri.mesh${task}" "${CMAKE_CURRENT_BINARY_DIR}/feti_beam_fine_tri.msh_${task}"
)
endforeach()


execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "/home/phuschke/nuto_project/nuto/myNutoExamples/meshFiles/2d/feti/feti_linear_benchmark.msh" "${CMAKE_CURRENT_BINARY_DIR}/linear_benchmark.msh"
)



foreach(file
        2d_gradient_three_point_bending
        2d_gradient_4_point_bending_symmetry
        2d_linear_feti_benchmark
        )

add_executable(${file} ${file}.cpp)
target_link_libraries(${file} ${MPI_LIBRARIES})
target_link_libraries(${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
target_link_libraries(${file} NuToVisualize)
target_link_libraries(${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endforeach()

