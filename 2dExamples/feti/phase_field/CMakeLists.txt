set(MESH_DIR "${CMAKE_SOURCE_DIR}/applications/custom/myNutoExamples/meshFiles/2d/feti/virtual_constraints")


if (ENABLE_MPI)
    foreach(task RANGE 0 1)
        foreach(mesh
                2d_single_edge_notched_tension_test_2_subdomains.mesh
                )

            execute_process(
                    COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${task}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${task}"
            )
        endforeach()
    endforeach()


    foreach(file
            2d_feti_structure_phase_field_single_edge_notched_tension_test
            )

        add_executable(${file} ${file}.cpp)
        target_link_libraries(${file} ${MPI_LIBRARIES})
        target_link_libraries(${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
        target_link_libraries(${file} NuToVisualize)
        target_link_libraries(${file} NuToFeti)
        target_link_libraries(${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endforeach()
endif ()
