add_subdirectory(meshes)
add_subdirectory(octaveScripts)


foreach(file
        2d_single_edge_notched_tension_local_damage
        )
    add_executable        (${file} ${file}.cpp)
    target_link_libraries (${file} Mechanics Math Base ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
    target_link_libraries (${file} Visualize)
    target_link_libraries (${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endforeach()