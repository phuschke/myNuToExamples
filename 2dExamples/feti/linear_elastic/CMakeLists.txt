set(MESH_DIR "${CMAKE_SOURCE_DIR}/applications/custom/myNutoExamples/meshFiles/2d/feti/virtual_constraints")

foreach(mesh
        feti_linear_benchmark.msh
        )
    execute_process(
            COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}"
    )
endforeach()

########################################################################################################################
########################################################################################################################

if (ENABLE_MPI)
    foreach(file
            2d_feti_structure
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
        2d_feti_structure_compare
        )
    add_executable(${file} ${file}.cpp)
    target_link_libraries(${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    target_link_libraries(${file} NuToVisualize)
    target_link_libraries(${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endforeach()