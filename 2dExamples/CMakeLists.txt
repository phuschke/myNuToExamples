#ADD_SUBDIRECTORY(interface)

ADD_SUBDIRECTORY(feti)
ADD_SUBDIRECTORY(phase_field)
ADD_SUBDIRECTORY(benchmark)
ADD_SUBDIRECTORY(gradient_damage)

set(MESH_DIR "${CMAKE_SOURCE_DIR}/applications/custom/myNutoExamples/meshFiles/2d")


execute_process(
  COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${MESH_DIR}/feti/fine/feti_beam.msh" "${CMAKE_CURRENT_BINARY_DIR}/feti_beam.msh"
)
