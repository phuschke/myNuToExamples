
set(MESH_DIR "/home/phuschke/nuto_project/nuto/applications/custom/myNutoExamples/meshFiles")

execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/3d/3d_3_point_bending.msh" "${CMAKE_CURRENT_BINARY_DIR}/3d_3_point_bending.msh")


foreach(file
        3d_3_point_bending
        3d_3_point_bending_test
        3d_3_point_bending_cuboid
        )

    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/3d/${file}.msh" "${CMAKE_CURRENT_BINARY_DIR}/${file}.msh")
endforeach()

foreach(file
        3d_linear_constraints
        3d_reference
        3d_singular_matrix
        )

    add_executable(${file} ${file}.cpp)
    target_link_libraries(${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    target_link_libraries(${file} NuToVisualize)
    target_link_libraries(${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endforeach()