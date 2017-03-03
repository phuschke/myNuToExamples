set(MESH_DIR "/home/phuschke/nuto_project/nuto/applications/custom/myNutoExamples/meshFiles/2d/feti/virtual_constraints")
foreach(task RANGE 0 3)
    foreach(mesh
            feti_beam_4_subdomains.mesh
            feti_beam_4_subdomains_60000_ele.mesh
            feti_beam_4_subdomains_1500000_ele.mesh
            )

        execute_process(
                COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${task}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${task}"
        )
    endforeach()
endforeach()


foreach(task RANGE 0 1)
    foreach(mesh
            feti_beam_coarse_2_subdomains.mesh
            feti_beam_coarse_2_subdomains_24_ele.mesh
            feti_beam_coarse_2_subdomains_6_ele.mesh
            )

        execute_process(
                COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${task}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${task}"
        )
    endforeach()
endforeach()


foreach(mesh
        feti_beam_coarse_2_subdomains.msh
        feti_beam_coarse_2_subdomains_24_ele_compare.msh
        feti_beam_coarse_2_subdomains_6_ele_compare.msh
         )

    execute_process(
            COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}"
    )
endforeach()


if (ENABLE_MPI)
foreach(file
        2d_feti_structure
        #2d_feti_structure_local_damage
        #2d_feti_structure_gradient_damage
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
        2d_feti_structure_compare
#        2d_feti_structure_local_damage_compare
#        2d_feti_structure_gradient_damage_compare
        )

    add_executable(${file} ${file}.cpp)
    target_link_libraries(${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    target_link_libraries(${file} NuToVisualize)
    target_link_libraries(${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endforeach()
