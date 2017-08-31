add_subdirectory(paraviewScripts)
add_subdirectory(octaveScripts)
#set(MESH_DIR "${CMAKE_SOURCE_DIR}/applications/custom/myNutoExamples/meshFiles/2d/feti/")
#
#########################################################################################################################
#########################################################################################################################
#
#foreach(task RANGE 0 5)
#    foreach(mesh
#            feti_4_point_bending_test_6_subs.mesh
#            )
#
#        execute_process(
#                COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${task}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${task}"
#        )
#    endforeach()
#endforeach()
#
#
#########################################################################################################################
#########################################################################################################################
#
#
#
#foreach(task RANGE 0 1)
#    foreach(mesh
#            feti_4_point_bending_test.mesh
#            )
#
#        execute_process(
#                COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}${task}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}${task}"
#        )
#    endforeach()
#endforeach()
#

if (ENABLE_MPI)
    foreach(file
            strong_scaling
            )

        add_executable(${file} ${file}.cpp)
        target_link_libraries(${file} ${MPI_LIBRARIES})
        target_link_libraries(${file} Mechanics Math Base ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
        target_link_libraries(${file} Visualize)
        target_link_libraries(${file} FetiLib)
        target_link_libraries(${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endforeach()
endif ()