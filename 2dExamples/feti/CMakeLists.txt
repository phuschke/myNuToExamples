


# create symlinks to the necessary mesh files
foreach(task RANGE 0 4)
    message("Task ${task}")
execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "/home/phuschke/gmshReader_project/build/examples/fourPointBending.mesh${task}" "${CMAKE_CURRENT_BINARY_DIR}/fourPointBending.msh_${task}"
)
  endforeach()

# 2d_feti_structure_phase_field_json
set(file 2d_gradient_4_point_bending_symmetry)
add_executable(${file} ${file}.cpp)
target_link_libraries(${file} ${MPI_LIBRARIES})
target_link_libraries(${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
target_link_libraries(${file} NuToVisualize)
target_link_libraries(${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

IF(MINGW)
  SET_TARGET_PROPERTIES(2d_feti_structure_phase_field_json PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
ENDIF(MINGW)


