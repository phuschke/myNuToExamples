set(MESH_DIR "${CMAKE_SOURCE_DIR}/applications/custom/myNutoExamples/meshFiles/2d/feti/virtual_constraints")

foreach(mesh
        feti_linear_benchmark.msh
        feti_linear_benchmark_2_subs_500k_nodes.msh
        feti_linear_benchmark_240k_nodes.msh
        )
    execute_process(
            COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}"
    )
endforeach()

########################################################################################################################
########################################################################################################################

foreach(task RANGE 0 1)
    foreach(mesh
            feti_beam_coarse_2_subdomains.mesh
            feti_beam_coarse_2_subdomains_24_ele.mesh
            feti_beam_coarse_2_subdomains_6_ele.mesh
            feti_linear_benchmark_2_subs.mesh
            feti_linear_benchmark_2_subs_500k_nodes.mesh
            feti_linear_benchmark_2_subs_240k_nodes.mesh
            )

        execute_process(
                COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${task}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${task}"
        )
    endforeach()
endforeach()

########################################################################################################################
########################################################################################################################

foreach(task RANGE 0 3)
    foreach(mesh
            feti_linear_benchmark_4_subs_500k_nodes.mesh
            feti_linear_benchmark_4_subs_240k_nodes.mesh
            )

        execute_process(
                COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${task}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${task}"
        )
    endforeach()
endforeach()

########################################################################################################################
########################################################################################################################

foreach(task RANGE 0 7)
    foreach(mesh
            feti_linear_benchmark_8_subs_500k_nodes.mesh
            feti_linear_benchmark_8_subs_240k_nodes.mesh
            )

        execute_process(
                COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${task}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${task}"
        )
    endforeach()
endforeach()

########################################################################################################################
########################################################################################################################

if (ENABLE_MPI)
    foreach(file
#            2d_feti_structure
            2d_feti_structure_manual_partitions
#            2d_feti_structure_projected_gmres
            )

        add_executable(${file} ${file}.cpp)
        target_link_libraries(${file} ${MPI_LIBRARIES})
        target_link_libraries(${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
        target_link_libraries(${file} NuToVisualize)
        target_link_libraries(${file} NuToFeti)
        target_link_libraries(${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endforeach()
endif ()

########################################################################################################################
########################################################################################################################

foreach(file
        #2d_feti_structure_compare
        )
    add_executable(${file} ${file}.cpp)
    target_link_libraries(${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    target_link_libraries(${file} NuToVisualize)
    target_link_libraries(${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endforeach()