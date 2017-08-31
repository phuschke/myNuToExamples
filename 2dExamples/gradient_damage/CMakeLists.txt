add_subdirectory(meshes)
add_subdirectory(octaveScripts)

#set(MESH_DIR "${CMAKE_SOURCE_DIR}/applications/custom/myNutoExamples/meshFiles/2d")
#
#
#foreach(mesh
#        2d_L_shaped_panel.msh
#        )
#  execute_process(
#          COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${mesh}" "${CMAKE_CURRENT_BINARY_DIR}/${mesh}"
#  )
#endforeach()
#
#
#
foreach(file
        2d_single_edge_notched_tension_test_quads_ged
        2d_test_convergence_of_damage_laws
        )
  add_executable        (${file} ${file}.cpp)
  target_link_libraries (${file} Mechanics Math Base ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
  target_link_libraries (${file} Visualize)
  target_link_libraries (${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endforeach()