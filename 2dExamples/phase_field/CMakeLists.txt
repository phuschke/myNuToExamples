add_subdirectory(meshes)
add_subdirectory(octaveScripts)
set(MESH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/meshes/")


foreach(file
        wedge_split_test.msh
        single_edge_notched_tension_test_quads.msh
        single_edge_notched_tension_test_quads_coarse.msh
        2d_single_edge_notched_tension_test_quads.msh
        )
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${file}" "${CMAKE_CURRENT_BINARY_DIR}/${file}")
endforeach()



foreach(file
        wedge_split_test
        single_edge_notched_tension_test_quads
        single_edge_notched_tension_test_quads_coarse
        )
    add_executable        (${file} ${file}.cpp)
    target_link_libraries (${file} Mechanics Math Base ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    target_link_libraries (${file} Visualize)
    target_link_libraries (${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endforeach()




