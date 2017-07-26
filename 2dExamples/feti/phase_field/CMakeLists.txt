add_subdirectory(paraviewScripts)
add_subdirectory(octaveScripts)
add_subdirectory(meshes)

set(MESH_DIR "${CMAKE_SOURCE_DIR}/applications/custom/myNutoExamples/meshFiles/2d/feti/phase_field")
if (ENABLE_MPI)
    foreach(task RANGE 0 1)
        foreach(mesh
                2d_single_edge_notched_tension_test_2_subdomains.mesh
                2d_single_edge_notched_tension_test_2_subdomains_quad.mesh
                )

            execute_process(
                    COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${task}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${task}"
            )
        endforeach()
    endforeach()


    foreach(task RANGE 0 3)
        foreach(mesh
                2d_single_edge_notched_tension_test_4_subdomains.mesh
                )

            execute_process(
                    COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${task}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${task}"
            )
        endforeach()
    endforeach()

    foreach(task RANGE 0 4)
        foreach(mesh
            2d_single_edge_notched_tension_test_quads_5_subdomains.mesh
            )

            execute_process(
                    COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${task}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${task}"
            )
        endforeach()
    endforeach()

    foreach(task RANGE 0 5)
        foreach(mesh
                2d_single_edge_notched_tension_test_6_subdomains.mesh
                )
            execute_process(
                    COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${task}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${task}"
            )
        endforeach()
    endforeach()

    foreach(task RANGE 0 2)
        foreach(mesh
                2d_single_edge_notched_tension_test_3_subdomains.mesh
                )

            execute_process(
                    COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${task}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${task}"
            )
        endforeach()
    endforeach()

    foreach(file
            2d_feti_structure_phase_field_single_edge_notched_tension_test
            feti_wedge_split_test
            )

        add_executable(${file} ${file}.cpp)
        target_link_libraries(${file} ${MPI_LIBRARIES})
        target_link_libraries(${file} Mechanics Math Base ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
        target_link_libraries(${file} Visualize)
        target_link_libraries(${file} FetiLib)
        target_link_libraries(${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endforeach()
endif ()
