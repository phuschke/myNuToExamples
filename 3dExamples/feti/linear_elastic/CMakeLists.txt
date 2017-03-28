set(MESH_DIR "${CMAKE_SOURCE_DIR}/applications/custom/myNutoExamples/meshFiles/3d/3d_feti/")

if(ENABLE_MPI)


    foreach(task RANGE 0 4)
        foreach(mesh
                3d_feti_5_subdomains_14688_ele.mesh
                3d_feti_5_subdomains_97768_ele.mesh
                3d_feti_5_subdomains_53328_ele.mesh
                )
            execute_process(
                    COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${task}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${task}"
            )
        endforeach()
    endforeach()




    foreach(file
            3d_feti_linear_elastic
            )

        add_executable(${file} ${file}.cpp)
        target_link_libraries(${file} ${MPI_LIBRARIES})
        target_link_libraries(${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
        target_link_libraries(${file} NuToVisualize)
        target_link_libraries(${file} NuToFeti)
        target_link_libraries(${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endforeach()

endif()