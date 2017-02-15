set(MESH_DIR "/home/phuschke/nuto_project/nuto/applications/custom/myNutoExamples/meshFiles")



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

execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/2d/feti/virtual_constraints/feti_beam_coarse_2_subdomains_24_ele.mesh${task}" "${CMAKE_CURRENT_BINARY_DIR}/feti_beam_coarse_2_subdomains_24_ele.mesh${task}"
)

execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/2d/feti/virtual_constraints/feti_beam_coarse_2_subdomains_6_ele.mesh${task}" "${CMAKE_CURRENT_BINARY_DIR}/feti_beam_coarse_2_subdomains_6_ele.mesh${task}"
)
endforeach()


execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/2d/feti/virtual_constraints/feti_beam_coarse_2_subdomains.msh" "${CMAKE_CURRENT_BINARY_DIR}/feti_beam_coarse_2_subdomains.msh"
)

execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/2d/feti/virtual_constraints/feti_beam_coarse_2_subdomains_compare_24_ele.msh" "${CMAKE_CURRENT_BINARY_DIR}/feti_beam_coarse_2_subdomains_compare_24_ele.msh"
)

execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/2d/feti/virtual_constraints/feti_beam_coarse_2_subdomains_compare_6_ele.msh" "${CMAKE_CURRENT_BINARY_DIR}/feti_beam_coarse_2_subdomains_compare_6_ele.msh"
)

if (ENABLE_MPI)
foreach(file
        #2d_feti_structure
        #2d_feti_structure_local_damage
        2d_feti_structure_gradient_damage
        )

    add_executable(${file} ${file}.cpp)
    target_link_libraries(${file} ${MPI_LIBRARIES})
    target_link_libraries(${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    target_link_libraries(${file} NuToVisualize)
    target_link_libraries(${file} NuToFeti)
    target_link_libraries(${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endforeach()

endif ()

foreach(file
#        2d_feti_structure_local_damage_compare
#        2d_feti_structure_gradient_damage_compare
        )

    add_executable(${file} ${file}.cpp)
    target_link_libraries(${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    target_link_libraries(${file} NuToVisualize)
    target_link_libraries(${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endforeach()
