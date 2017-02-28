
set(MESH_DIR "${CMAKE_SOURCE_DIR}/applications/custom/myNutoExamples/meshFiles")
execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/2d/2d_single_edge_notched_tension_test.msh" "${CMAKE_CURRENT_BINARY_DIR}/2d_single_edge_notched_tension_test.msh")


foreach(file
        2d_phase_field_single_edge_notched_tension
        2d_phase_field_3_point_miehe
        )
    add_executable        (${file} ${file}.cpp)
    target_link_libraries (${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    target_link_libraries (${file} NuToVisualize)
    target_link_libraries (${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endforeach()




