set(MESH_DIR "/home/phuschke/nuto_project/nuto/myNutoExamples/meshFiles")



foreach(task RANGE 0 3)

    execute_process(
        COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/2d/feti/virtual_constraints/feti_beam_4_subdomains.mesh${task}" "${CMAKE_CURRENT_BINARY_DIR}/feti_beam_4_subdomains.mesh${task}"
        )

    execute_process(
        COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/2d/feti/virtual_constraints/feti_beam_4_subdomains_60000_ele.mesh${task}" "${CMAKE_CURRENT_BINARY_DIR}/feti_beam_4_subdomains_60000_ele.mesh${task}"
        )

    execute_process(
        COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/2d/feti/virtual_constraints/feti_beam_4_subdomains_1500000_ele.mesh${task}" "${CMAKE_CURRENT_BINARY_DIR}/feti_beam_4_subdomains_1500000_ele.mesh${task}"
        )

endforeach()


foreach(task RANGE 0 1)
    message("Task ${task}")
execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/2d/feti/virtual_constraints/feti_beam_coarse_2_subdomains.mesh${task}" "${CMAKE_CURRENT_BINARY_DIR}/feti_beam_coarse_2_subdomains.mesh${task}"
)
endforeach()





foreach(file
        2d_feti_structure
        )

    add_executable(${file} ${file}.cpp)
    target_link_libraries(${file} ${MPI_LIBRARIES})
    target_link_libraries(${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    target_link_libraries(${file} NuToVisualize)
    target_link_libraries(${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endforeach()
