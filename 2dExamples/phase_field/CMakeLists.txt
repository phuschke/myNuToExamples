#
#set(MESH_DIR "${CMAKE_SOURCE_DIR}/applications/custom/myNutoExamples/meshFiles/2d")
#
#
#foreach(file
#2d_single_edge_notched_tension_test.msh
#2d_single_edge_notched_tension_test_quads.msh
#2d_miehe_symmetric_three_point_bending.msh
#2d_miehe_symmetric_three_point_bending_no_notch.msh
#2d_matrix_thought_experiment.msh
#2d_L_shaped_panel.msh
#
#)
#execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/${file}" "${CMAKE_CURRENT_BINARY_DIR}/${file}")
#endforeach()
#
#
#
#foreach(file
#        2d_phase_field_single_edge_notched_tension
#        #2d_phase_field_3_point_miehe
#        #2d_phase_field_tensile_test
#        #2d_phase_field_L_shaped_panel
#        )
#    add_executable        (${file} ${file}.cpp)
#    target_link_libraries (${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
#    target_link_libraries (${file} NuToVisualize)
#    target_link_libraries (${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
#endforeach()
#
#
#
#
