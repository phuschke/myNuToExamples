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
#foreach(file
#        2d_gradient_L_shaped_panel
#        )
#  add_executable        (${file} ${file}.cpp)
#  target_link_libraries (${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
#  target_link_libraries (${file} NuToVisualize)
#  target_link_libraries (${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
#endforeach()